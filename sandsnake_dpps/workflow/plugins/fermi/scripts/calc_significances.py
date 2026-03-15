import numpy as np
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.table import QTable
from astropy.time import Time
import astropy.units as u
from pathlib import Path
import matplotlib
import matplotlib.pyplot as plt
import sys
from scipy.optimize import curve_fit
from gammapy.irf import load_irf_dict_from_file
from gammapy.data import Observation, Observations
from gammapy.maps import MapAxis, RegionGeom
from gammapy.datasets import SpectrumDataset, SpectrumDatasetOnOff, Datasets
from gammapy.makers import SpectrumDatasetMaker, SafeMaskMaker
from gammapy.modeling.models import (
    PowerLawSpectralModel,
    LogParabolaSpectralModel,
    SuperExpCutoffPowerLaw4FGLDR3SpectralModel,
    ExpCutoffPowerLawNormSpectralModel,
    EBLAbsorptionNormSpectralModel,
    SkyModel,
)

from collections import defaultdict
from regions import CircleSkyRegion
import re
import logging

if matplotlib.get_backend() == "pgf":
    from matplotlib.backends.backend_pgf import PdfPages
else:
    from matplotlib.backends.backend_pdf import PdfPages

from core.scripts.mc.irf_plots import add_sensitivity_comparisons
from .process_catalog import (
    calc_delta_B,
    prod_site_B_declination,
    prod_site_B_inclination,
)

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

obstimes = [0.5, 5.0, 50.0]
ebl_model = "dominguez"
ctao_north = EarthLocation.of_site("Roque de los Muchachos")
offset = 0.7 * u.deg
on_radius = 0.35 * u.deg  # doesnt matter for pointlike IRFs with cuts already applied
n_off_regions = 5
n_obs = 100


def parse_irf_path(path: str):
    m = re.search(r"zen_(\d+)/az_(\d+)/.*obs_([0-9]+(?:\.[0-9]+)?)_hours", path)
    if m is None:
        raise ValueError(path)

    zen = int(m.group(1))
    az = int(m.group(2))
    obstime = m.group(3)

    return zen, az, obstime


def create_dict(paths):
    dd = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))

    for p in paths:
        zen, az, obstime = parse_irf_path(p)
        dd[zen][az][obstime].append(p)

    return dd


def get_zen_az_pairs(zen_az_dict):
    return sorted((zen, az) for zen, az_dict in zen_az_dict.items() for az in az_dict)


def delta_B_mc_grid(pointings):
    frame = AltAz(location=ctao_north, obstime=Time.now())
    table = QTable(names=("zen", "az", "delta_B"), units=(u.deg, u.deg, u.deg))

    for zen, az in pointings:
        alt = 90 * u.deg - zen
        pointing = SkyCoord(alt=alt, az=az, frame=frame)
        delta_B = calc_delta_B(
            pointing, None, prod_site_B_declination, prod_site_B_inclination
        )

        table.add_row((zen, az, delta_B))

    return table


def get_nearest_node(grid_table, alt_median, delta_B_median):
    zen_target = 90 * u.deg - alt_median

    zeniths = np.unique(grid_table["zen"])
    nearest_zen = zeniths[np.argmin(np.abs(zeniths - zen_target))]

    cand = grid_table[grid_table["zen"] == nearest_zen]
    idx = np.argmin(np.abs(cand["delta_B"] - delta_B_median))
    return (cand[idx]["zen"], cand[idx]["az"])


def load_irfs(irfs_dict, bench_dict, nn_pointing):
    zen, az = nn_pointing

    nn_irfs_dict = defaultdict()
    for obstime, path in irfs_dict[zen][az].items():
        nn_irfs_dict[obstime] = load_irf_dict_from_file(path)

    nn_sens_dict = defaultdict()
    for obstime, path in bench_dict[zen][az].items():
        nn_sens_dict[obstime] = QTable.read(path, hdu="SENSITIVITY")

    return nn_irfs_dict, nn_sens_dict


def get_irfs_sens(irfs, benchmarks, source_row):
    irfs_dict = create_dict(irfs)
    bench_dict = create_dict(benchmarks)
    mc_pointings = get_zen_az_pairs(irfs_dict)
    delta_B_table = delta_B_mc_grid(mc_pointings)
    nn_pointing = get_nearest_node(
        delta_B_table, source_row["alt_median"], source_row["delta_B_median"]
    )
    nn_irfs_dict, nn_sens_dict = load_irfs(irfs_dict, bench_dict, nn_pointing)

    return nn_irfs_dict, nn_sens_dict


def create_spectral_model(source_table):
    spec_type = source_table["SpectrumType"]
    sed_class = None
    if "SED_class" in source_table.colnames:
        sed_class = source_table["SED_class"]

    # TODO: Implement Redshift bounds for sources without redshift -> 0 - 0.5 ?
    redshift = 0
    if "Redshift" in source_table.colnames:
        if np.isfinite(source_table["Redshift"]):
            redshift = source_table["Redshift"]

    # Spectral Model
    if "FHL" in source_table["Source_Name"]:
        if spec_type == "PowerLaw":
            spec_model = PowerLawSpectralModel(
                amplitude=source_table["Flux_Density"] / u.ph,
                reference=source_table["Pivot_Energy"],
                index=source_table["PowerLaw_Index"],
            )
        elif spec_type == "LogParabola":
            spec_model = LogParabolaSpectralModel(
                amplitude=source_table["Flux_Density"] / u.ph,
                reference=source_table["Pivot_Energy"],
                alpha=source_table["Spectral_Index"],
                beta=source_table["beta"],
            )
        else:
            log.critical(f"Spectral Model: {spec_type} not implemented!")
            sys.exit(1)

    else:
        if spec_type == "PowerLaw":
            spec_model = PowerLawSpectralModel(
                amplitude=source_table["PL_Flux_Density"] / u.ph,
                reference=source_table["Pivot_Energy"],
                index=source_table["PL_Index"],
            )
        elif spec_type == "LogParabola":
            spec_model = LogParabolaSpectralModel(
                amplitude=source_table["LP_Flux_Density"] / u.ph,
                reference=source_table["Pivot_Energy"],
                alpha=source_table["LP_Index"],
                beta=source_table["LP_beta"],
            )

        # TODO: Better force LogParabola instead?
        elif spec_type == "PLSuperExpCutoff":
            spec_model = SuperExpCutoffPowerLaw4FGLDR3SpectralModel(
                amplitude=source_table["PLEC_Flux_Density"] / u.ph,
                reference=source_table["Pivot_Energy"],
                index_1=source_table["PLEC_IndexS"],
                index_2=source_table["PLEC_Exp_Index"],
                expfactor=source_table["PLEC_ExpfactorS"],
            )
        else:
            log.critical(f"Spectral Model: {spec_type} not implemented!")
            sys.exit(1)

    # TODO: Cutoff überhaupt sinvoll, wenn kein AGN? Und PL & LP gut gefittet?
    # Energy Cutoff
    cutoff_spec_model = None
    if spec_type in ["PowerLaw", "LogParabola"]:
        # or 10 TeV for all?
        if sed_class in ["LSP", "ISP"]:
            cutoff_energy = 0.1 * u.TeV
        elif sed_class == "HSP":
            cutoff_energy = 1 * u.TeV
        else:
            cutoff_energy = 10 * u.TeV

        # corrections to observers frame
        cutoff_energy = cutoff_energy / (1 + redshift)

        cutoff_spec_model = ExpCutoffPowerLawNormSpectralModel(
            norm=1,
            index=0,
            lambda_=1 / cutoff_energy,
            alpha=1,
            reference=cutoff_energy / 10,
        )

    # EBL Absorption
    ebl_abs_model = EBLAbsorptionNormSpectralModel.read_builtin(
        ebl_model, redshift=redshift
    )

    model = spec_model * ebl_abs_model
    if cutoff_spec_model is not None:
        model = model * cutoff_spec_model

    return SkyModel(spectral_model=model, name=source_table["Source_Name"])


def index_observations_by_obstime(observations: Observations):
    idx = defaultdict(list)
    for obs in observations:
        idx[obs.meta["obstime"]].append(obs)
    return {k: Observations(v) for k, v in idx.items()}


def get_observations(irfs_dict, pointing):
    obs_list = []

    for obs_id, (obstime, irfs) in enumerate(irfs_dict.items(), start=1):
        obs = Observation.create(
            obs_id=obs_id,
            pointing=pointing,
            livetime=float(obstime) * u.h,
            irfs=irfs,
            location=ctao_north,
        )
        obs.meta["obstime"] = obstime  # str
        obs_list.append(obs)

    return index_observations_by_obstime(Observations(obs_list))


def create_spectrum_dataset_onoff(observation, on_region, source_model):
    energy_axis_reco = observation.bkg.axes["energy"]
    energy_axis_true = MapAxis.from_energy_bounds(
        0.3 * energy_axis_reco.edges[0],
        3.0 * energy_axis_reco.edges[-1],
        nbin=3 * len(energy_axis_reco.edges),
        name="energy_true",
    )

    geom = RegionGeom.create(region=on_region, axes=[energy_axis_reco])
    dataset_empty = SpectrumDataset.create(
        geom=geom, energy_axis_true=energy_axis_true, name="obs"
    )
    maker = SpectrumDatasetMaker(
        containment_correction=False,
        use_region_center=True,
        selection=["exposure", "edisp", "background"],
    )
    safe_mask_maker = SafeMaskMaker(methods=["bkg-peak"])

    # TODO: Create galactic lat dependent scaling factor for background?
    # dataset.background.data *= scaling_factor
    dataset = maker.run(dataset_empty, observation)
    dataset = safe_mask_maker.run(dataset, observation)

    dataset.models = source_model

    dataset_on_off = SpectrumDatasetOnOff.from_spectrum_dataset(
        dataset=dataset, acceptance=1, acceptance_off=n_off_regions
    )

    return dataset_on_off


def fake_data(dataset_on_off, source_model):
    datasets = Datasets()
    for idx in range(n_obs):
        ds = dataset_on_off.copy(name=f"obs_{idx}")
        # Dataset.copy() doesnt really work with its models - need to copy model by hand
        ds.models = source_model.copy()
        ds.fake(random_state=idx, npred_background=ds.npred_background())
        ds.meta_table["OBS_ID"] = [idx]
        datasets.append(ds)

    table = datasets.info_table()
    sigma = table["sqrt_ts"]

    return table, sigma


def plot_source_model(sens_dict, source_model, out_dir):
    # TODO: Add PROD5 Sens?
    fig, ax = plt.subplots()
    e_lim = [5.0e-3, 5.0e2]

    for obstime in obstimes:
        sens = sens_dict[str(obstime)]

        ax.errorbar(
            (0.5 * (sens["ENERG_LO"] + sens["ENERG_HI"])).flatten(),
            sens["ENERGY_FLUX_SENSITIVITY"].flatten(),
            xerr=0.5 * (sens["ENERG_HI"] - sens["ENERG_LO"]),
            ls="",
            label=f"CTAO-N - {obstime}h",
        )

    # Plot Source Model
    ax = source_model.spectral_model.plot(
        energy_bounds=[0.01, 100] * u.TeV,
        ax=ax,
        label=source_model.name,
        sed_type="e2dnde",
    )

    add_sensitivity_comparisons(ax, energy_limits=e_lim)

    ax.set_ylim(3.0e-14, 1.0e-9)
    ax.set_xlim(e_lim)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_title("Sensitivity")
    ax.set_xlabel(r"$E_{True}$ / TeV")
    ax.set_ylabel(r"$E^{2} \times$ Flux Sensitivity / $erg \cdot cm^{-2} \cdot s^{1}$")
    ax.grid(which="both", linestyle=":")
    ax.legend(loc="upper right", fontsize="small")

    fig.tight_layout()

    with PdfPages(out_dir / f"{source_model.name}_sensitivity.pdf") as pdf:
        pdf.savefig(fig)
        plt.close()


def predict_obstime(sigma_dict, sigma_target, irfs_dict):
    obstime_str = sorted(irfs_dict.keys(), key=float)
    obstimes_array = np.array([float(o) for o in obstime_str], dtype=float)

    sigma_mean = np.array(
        [sigma_dict[f"{o}h"].mean() for o in obstime_str], dtype=float
    )
    sigma_std = np.array([sigma_dict[f"{o}h"].std() for o in obstime_str], dtype=float)

    def model(T, a):
        # S ~ sqrt(T)
        return a * np.sqrt(T)

    # Fit
    popt, pcov = curve_fit(
        model, obstimes_array, sigma_mean, sigma=sigma_std, absolute_sigma=True
    )
    a = popt[0]
    da = np.sqrt(pcov[0, 0])

    # Error prop: T = (S/a)^2  →  dT/T = 2 * da/a  (S fix)
    obstime_pred = (sigma_target / a) ** 2
    obstime_std = obstime_pred * 2.0 * (da / a)

    return obstime_pred, obstime_std


def main(source, output, irfs, benchmarks):
    source_table = QTable.read(source, format="ascii.ecsv")
    source_row = source_table[0]

    # TODO: interpolate IRFs?
    irfs_dict, sens_dict = get_irfs_sens(irfs, benchmarks, source_row)

    source_model = create_spectral_model(source_row)
    source_position = SkyCoord(
        source_row["RAJ2000"], source_row["DEJ2000"], unit="deg", frame="icrs"
    )
    pointing = SkyCoord(source_position.ra + offset, source_position.dec)
    on_region = CircleSkyRegion(center=source_position, radius=on_radius)

    fake_table_dict = defaultdict()
    fake_sigma_dict = defaultdict()
    # TODO:Add gamma diffusive background model for galactic sources?
    observations = get_observations(irfs_dict, pointing)

    spec_dataset_on_off_50h = create_spectrum_dataset_onoff(
        observations["50"], on_region, source_model
    )
    table_50h, sigma_50h = fake_data(spec_dataset_on_off_50h, source_model)
    fake_table_dict["50h"] = table_50h
    fake_sigma_dict["50h"] = sigma_50h

    if sigma_50h.mean() > 5:
        for obstime in [x for x in irfs_dict.keys() if x != "50"]:
            spec_dataset_on_off = create_spectrum_dataset_onoff(
                observations[f"{int(obstime)}"],
                on_region,
                source_model,
            )

            table, sigma = fake_data(spec_dataset_on_off, source_model)
            fake_table_dict[f"{obstime}h"] = table
            fake_sigma_dict[f"{obstime}h"] = sigma

        obstime_5sigma, obstime_5sigma_std = predict_obstime(
            fake_sigma_dict, 5.0, irfs_dict
        )

        plot_source_model(sens_dict, source_model, Path(output).parent)

        source_table["obstime_5s"] = obstime_5sigma
        source_table["obstime_5s_std"] = obstime_5sigma_std
        for obstime in irfs_dict.keys():
            source_table[f"sigma_{obstime}h_mean"] = fake_sigma_dict[
                f"{obstime}h"
            ].mean()
            source_table[f"sigma_{obstime}h_std"] = fake_sigma_dict[f"{obstime}h"].std()

    else:
        source_table["obstime_5s"] = np.nan
        source_table["obstime_5s_std"] = np.nan
        source_table["sigma_50.0h_mean"] = sigma_50h.mean()
        source_table["sigma_50.0h_std"] = sigma_50h.std()
        for obstime in ["0.5", "5.0"]:
            source_table[f"sigma_{obstime}h_mean"] = np.nan
            source_table[f"sigma_{obstime}h_std"] = np.nan

    # Write table to new output
    source_table.write(Path(output), format="ascii.ecsv", overwrite=True)


if __name__ == "__main__":
    smk = globals().get("snakemake")
    if smk is not None:
        main(
            smk.input.source,
            smk.output[0],
            smk.input.irfs,
            smk.input.benchmarks,
        )
    else:
        import argparse

        parser = argparse.ArgumentParser()
        parser.add_argument("--source", required=True)
        parser.add_argument("--output", required=True)
        parser.add_argument("--irfs", nargs="+", required=True)
        parser.add_argument("--benchmarks", nargs="+", required=True)
        args = parser.parse_args()

        main(args.source, args.output, args.irfs, args.benchmarks)
