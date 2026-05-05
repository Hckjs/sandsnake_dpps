from __future__ import annotations

import numpy as np
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.table import QTable
from astropy.time import Time
import astropy.units as u
from pathlib import Path
import matplotlib.pyplot as plt
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
from dataclasses import dataclass
from typing import Iterable, Any
from regions import CircleSkyRegion
import re

from core.scripts.mc.irf_plots import add_sensitivity_comparisons
from .process_catalog import (
    VisibilityConfig,
    get_B_direction,
)

EXTRAGALACTIC_CLASSES = {
    "agn",
    "bcu",
    "bll",
    "fsrq",
    "rdg",
    "sey",
    "nlsy1",
    "ssrq",
    "sbg",
    "css",
}

GALACTIC_CLASSES = {
    "psr",
    "msp",
    "pwn",
    "snr",
    "spp",
    "glc",
    "hmb",
    "lmb",
    "bin",
    "nov",
    "sfr",
    "mc",
    "hii",
}

UNCERTAIN_CLASSES = {
    "",
    "--",
    "nan",
    "none",
    "unk",
    "unid",
    "unassociated",
}


@dataclass(slots=True)
class IRFNode:
    site_name: str
    location: EarthLocation
    prod_site_B_declination: u.Quantity
    prod_site_B_inclination: u.Quantity
    frame: AltAz

    zen: u.Quantity
    az: u.Quantity
    alt: u.Quantity
    pointing: SkyCoord
    delta_b: u.Quantity
    sin_delta: float
    cos_theta: float

    obstime_paths_irfs: dict[float, Path]
    obstime_paths_benchmarks: dict[float, Path]

    @property
    def key(self) -> tuple[int, int]:
        return (
            int(self.zen.to_value(u.deg)),
            int(self.az.to_value(u.deg)),
        )

    @property
    def obstimes_irfs(self) -> list[float]:
        return sorted(self.obstime_paths_irfs)

    @property
    def obstimes_benchmarks(self) -> list[float]:
        return sorted(self.obstime_paths_benchmarks)

    @property
    def common_obstimes(self) -> list[float]:
        return sorted(
            set(self.obstime_paths_irfs).intersection(self.obstime_paths_benchmarks)
        )

    def get_irf_path(self, obstime: float) -> Path:
        return self.obstime_paths_irfs[float(obstime)]

    def get_benchmark_path(self, obstime: float) -> Path:
        return self.obstime_paths_benchmarks[float(obstime)]

    def load_benchmark(self, obstime: float) -> QTable:
        path = self.get_benchmark_path(obstime)
        return QTable.read(path, hdu="SENSITIVITY")

    def load_irfs(self, obstime: float):
        path = self.get_irf_path(obstime)
        return load_irf_dict_from_file(path)

    def load_observations(self, start_obs_id: int = 1) -> Observations:
        observations = []
        for obs_id, obstime in enumerate(self.common_obstimes, start=start_obs_id):
            obs = Observation.create(
                obs_id=obs_id,
                pointing=self.pointing,
                livetime=float(obstime) * u.h,
                irfs=self.load_irfs(obstime),
                location=self.location,
            )
            obs.meta["obstime"] = str(obstime)
            observations.append(obs)

        return Observations(observations)


class SubarrayIRFs:
    """
    Collection of IRF and benchmark paths grouped into fixed (zen, az) nodes.

    The site-specific properties are initialized once and copied into every
    produced ``IRFNode`` dataclass instance.
    """

    # check common/paths.smk PATH["core:template:irfs"]
    # currently "/zen_{zen}/az_{az}/irfs_zen_{zen}_az_{az}_obs_{obstime}_hours.fits.gz"
    PATH_PATTERN = re.compile(
        r"/zen_(?P<zen>\d+)/az_(?P<az>\d+)/.*obs_(?P<obstime>[0-9]+(?:\.[0-9]+)?)_hours"
    )

    def __init__(
        self,
        irf_paths: Iterable[str | Path] | None = None,
        benchmark_paths: Iterable[str | Path] | None = None,
        *,
        site_name: str = "Roque de los Muchachos",
        b_declination: u.Quantity = VisibilityConfig.prod_site_B_declination,
        b_inclination: u.Quantity = VisibilityConfig.prod_site_B_inclination,
        frame_obstime: Time | None = None,
    ):
        self.site_name = site_name
        self.location = EarthLocation.of_site(site_name)
        self.prod_site_B_declination = b_declination
        self.prod_site_B_inclination = b_inclination
        self.frame_obstime = Time.now() if frame_obstime is None else frame_obstime
        self.frame = AltAz(location=self.location, obstime=self.frame_obstime)
        self.nodes: dict[tuple[int, int], IRFNode] = {}

        if irf_paths is not None or benchmark_paths is not None:
            self.nodes = self.build_nodes(
                irf_paths=() if irf_paths is None else irf_paths,
                benchmark_paths=() if benchmark_paths is None else benchmark_paths,
            )

    @classmethod
    def parse_input_path(
        cls, path: str | Path
    ) -> tuple[u.Quantity, u.Quantity, float, Path]:
        path = Path(path)
        match = cls.PATH_PATTERN.search(path.as_posix())
        if match is None:
            raise ValueError(f"Could not parse path: {path}")

        zen = float(match.group("zen")) * u.deg
        az = float(match.group("az")) * u.deg
        obstime = float(match.group("obstime"))
        return zen, az, obstime, path

    def _build_pointing(
        self, zen: u.Quantity, az: u.Quantity
    ) -> tuple[u.Quantity, SkyCoord]:
        alt = 90.0 * u.deg - zen
        pointing = SkyCoord(alt=alt, az=az, frame=self.frame)
        return alt, pointing

    def _calc_delta_b(self, pointing: SkyCoord) -> u.Quantity:
        b_vec = get_B_direction(
            self.prod_site_B_declination,
            self.prod_site_B_inclination,
        )
        los_vec = np.asarray(pointing.cartesian.xyz.value, dtype=float)
        cos_angle = np.clip(np.dot(los_vec, b_vec), -1.0, 1.0)
        angle = np.arccos(cos_angle) * u.rad
        return np.abs(90.0 * u.deg - angle.to(u.deg))

    def _group_paths_by_node_and_obstime(
        self,
        paths: Iterable[str | Path],
    ) -> dict[tuple[int, int], dict[float, Path]]:
        grouped: dict[tuple[int, int], dict[float, Path]] = defaultdict(dict)

        for path in paths:
            zen, az, obstime, parsed_path = self.parse_input_path(path)
            key = (int(zen.to_value(u.deg)), int(az.to_value(u.deg)))

            if obstime in grouped[key]:
                raise ValueError(
                    f"Duplicate obstime={obstime} for node {key}: "
                    f"{grouped[key][obstime]} and {parsed_path}"
                )

            grouped[key][obstime] = parsed_path

        return dict(grouped)

    def create_node(
        self,
        irf_paths: dict[float, Path],
        benchmark_paths: dict[float, Path],
    ) -> IRFNode:
        if not irf_paths and not benchmark_paths:
            raise ValueError("Cannot build an IRFNode without any paths")

        ref_paths = irf_paths if irf_paths else benchmark_paths
        ref_obstime = next(iter(ref_paths))
        zen_ref, az_ref, _, _ = self.parse_input_path(ref_paths[ref_obstime])

        for obstime, path in {**irf_paths, **benchmark_paths}.items():
            zen, az, _, _ = self.parse_input_path(path)
            if not u.isclose(zen, zen_ref) or not u.isclose(az, az_ref):
                raise ValueError(
                    "All paths passed to create_node() must belong to the same "
                    f"(zen, az) node, got {(zen_ref, az_ref)} and {(zen, az)}"
                )

        alt, pointing = self._build_pointing(zen_ref, az_ref)
        delta_b = self._calc_delta_b(pointing)
        sin_delta = float(np.sin(delta_b.to_value(u.rad)))
        cos_theta = float(np.cos(zen_ref.to_value(u.rad)))

        return IRFNode(
            site_name=self.site_name,
            location=self.location,
            prod_site_B_declination=self.prod_site_B_declination,
            prod_site_B_inclination=self.prod_site_B_inclination,
            frame=self.frame,
            zen=zen_ref,
            az=az_ref,
            alt=alt,
            pointing=pointing,
            delta_b=delta_b,
            sin_delta=sin_delta,
            cos_theta=cos_theta,
            obstime_paths_irfs=dict(sorted(irf_paths.items())),
            obstime_paths_benchmarks=dict(sorted(benchmark_paths.items())),
        )

    def build_nodes(
        self,
        irf_paths: Iterable[str | Path],
        benchmark_paths: Iterable[str | Path],
    ) -> dict[tuple[int, int], IRFNode]:
        grouped_irfs = self._group_paths_by_node_and_obstime(irf_paths)
        grouped_benchmarks = self._group_paths_by_node_and_obstime(benchmark_paths)

        all_keys = sorted(set(grouped_irfs) | set(grouped_benchmarks))

        return {
            key: self.create_node(
                irf_paths=grouped_irfs.get(key, {}),
                benchmark_paths=grouped_benchmarks.get(key, {}),
            )
            for key in all_keys
        }

    def get_nearest_node_key(
        self,
        cos_theta_mean: float,
        sin_delta_mean: float,
        weight_cos_theta: float = 1.0,
        weight_sin_delta: float = 1.0,
    ) -> tuple[int, int]:
        if not self.nodes:
            raise ValueError("No nodes available in SubarrayIRFs")

        def distance2(node: IRFNode) -> float:
            d_cos = node.cos_theta - cos_theta_mean
            d_sin = node.sin_delta - sin_delta_mean
            return weight_cos_theta * d_cos**2 + weight_sin_delta * d_sin**2

        best_node = min(self.nodes.values(), key=distance2)
        return best_node.key

    def get_nearest_node(
        self,
        cos_theta_mean: float,
        sin_delta_mean: float,
        weight_cos_theta: float = 1.0,
        weight_sin_delta: float = 1.0,
    ) -> IRFNode:
        return self.nodes[
            self.get_nearest_node_key(
                cos_theta_mean=cos_theta_mean,
                sin_delta_mean=sin_delta_mean,
                weight_cos_theta=weight_cos_theta,
                weight_sin_delta=weight_sin_delta,
            )
        ]


class Source:
    """
    Single source loaded from a one-row ECSV file.

    This class only contains intrinsic / catalog-derived source properties and the
    corresponding spectral model. No observation- or analysis-specific state is
    stored here.
    """

    def __init__(
        self,
        path: str | Path,
        *,
        ebl_model: str = "dominguez",
    ):
        self.path = Path(path)
        self.table = QTable.read(self.path, format="ascii.ecsv")

        if len(self.table) != 1:
            raise ValueError(
                f"Expected exactly one row in source table {self.path}, "
                f"got {len(self.table)}"
            )

        self.ebl_model = ebl_model

        self._row = self.table[0]
        self._position = SkyCoord(
            self._row["RAJ2000"],
            self._row["DEJ2000"],
            unit="deg",
            frame="icrs",
        )
        self._spectral_model = self._create_spectral_model(self._row)

    @property
    def row(self):
        return self._row

    @property
    def name(self) -> str:
        return str(self.row["Source_Name"])

    @property
    def position(self) -> SkyCoord:
        return self._position

    @property
    def spectral_model(self) -> SkyModel:
        return self._spectral_model

    def _create_spectral_model(self, row) -> SkyModel:
        spec_type = row["SpectrumType"]
        sed_class = row["SED_class"] if "SED_class" in row.colnames else None

        redshift = 0.0
        if "Redshift" in row.colnames and np.isfinite(row["Redshift"]):
            redshift = float(row["Redshift"])

        if "FHL" in str(row["Source_Name"]):
            if spec_type == "PowerLaw":
                spec_model = PowerLawSpectralModel(
                    amplitude=row["Flux_Density"] / u.ph,
                    reference=row["Pivot_Energy"],
                    index=row["PowerLaw_Index"],
                )
            elif spec_type == "LogParabola":
                spec_model = LogParabolaSpectralModel(
                    amplitude=row["Flux_Density"] / u.ph,
                    reference=row["Pivot_Energy"],
                    alpha=row["Spectral_Index"],
                    beta=row["beta"],
                )
            else:
                raise ValueError(
                    f"Spectral model {spec_type!r} not implemented for FHL source"
                )
        else:
            if spec_type == "PowerLaw":
                spec_model = PowerLawSpectralModel(
                    amplitude=row["PL_Flux_Density"] / u.ph,
                    reference=row["Pivot_Energy"],
                    index=row["PL_Index"],
                )
            elif spec_type == "LogParabola":
                spec_model = LogParabolaSpectralModel(
                    amplitude=row["LP_Flux_Density"] / u.ph,
                    reference=row["Pivot_Energy"],
                    alpha=row["LP_Index"],
                    beta=row["LP_beta"],
                )
            elif spec_type == "PLSuperExpCutoff":
                spec_model = SuperExpCutoffPowerLaw4FGLDR3SpectralModel(
                    amplitude=row["PLEC_Flux_Density"] / u.ph,
                    reference=row["Pivot_Energy"],
                    index_1=row["PLEC_IndexS"],
                    index_2=row["PLEC_Exp_Index"],
                    expfactor=row["PLEC_ExpfactorS"],
                )
            else:
                raise ValueError(f"Spectral model {spec_type!r} not implemented")

        cutoff_spec_model = None
        if spec_type in ["PowerLaw", "LogParabola"]:
            if sed_class in ["LSP", "ISP"]:
                cutoff_energy = 0.1 * u.TeV
            elif sed_class == "HSP":
                cutoff_energy = 1.0 * u.TeV
            else:
                cutoff_energy = 10.0 * u.TeV

            cutoff_energy = cutoff_energy / (1.0 + redshift)

            cutoff_spec_model = ExpCutoffPowerLawNormSpectralModel(
                norm=1.0,
                index=0.0,
                lambda_=1.0 / cutoff_energy,
                alpha=1.0,
                reference=cutoff_energy / 10.0,
            )

        ebl_abs_model = EBLAbsorptionNormSpectralModel.read_builtin(
            self.ebl_model,
            redshift=redshift,
        )

        model = spec_model * ebl_abs_model
        if cutoff_spec_model is not None:
            model = model * cutoff_spec_model

        return SkyModel(spectral_model=model, name=self.name)

    def append_columns(self, values: dict[str, Any]) -> None:
        for key, value in values.items():
            self.table[key] = [value]

    def write(self, output_path: str | Path, *, overwrite: bool = True) -> None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        self.table.write(output_path, format="ascii.ecsv", overwrite=overwrite)


class SourceAnalysis:
    """
    Concrete analysis of one source for one IRF node.

    This class combines:
    - a Source object
    - one IRFNode
    - observation-specific analysis settings
    """

    def __init__(
        self,
        source: Source,
        irf_node,
        *,
        offset: u.Quantity = 0.7 * u.deg,
        on_radius: u.Quantity = 0.35 * u.deg,
        n_off_regions: int = 5,
        n_fake_datasets: int = 100,
    ):
        self.source = source
        self.irf_node = irf_node

        self.offset = offset
        self.on_radius = on_radius
        self.n_off_regions = n_off_regions
        self.n_fake_datasets = n_fake_datasets

        self.pointing = SkyCoord(
            self.source.position.ra + self.offset,
            self.source.position.dec,
            frame="icrs",
        )
        self.on_region = CircleSkyRegion(
            center=self.source.position,
            radius=self.on_radius,
        )

        self.observations: Observations | None = None
        self.datasets_onoff: dict[float, SpectrumDatasetOnOff] | None = None
        self.fake_results: dict[float, dict[str, Any]] | None = None
        self.obstime_5s: float = np.nan
        self.obstime_5s_std: float = np.nan

    @staticmethod
    def _format_obstime(obstime: float) -> str:
        return f"{float(obstime):g}"

    def create_observations(
        self,
        *,
        obstimes: list[float] | None = None,
        start_obs_id: int = 1,
    ) -> Observations:
        if obstimes is None:
            obstimes = self.irf_node.common_obstimes

        if not obstimes:
            raise ValueError(
                f"No common observation times available for IRF node {self.irf_node.key}"
            )

        observations = []
        for obs_id, obstime in enumerate(sorted(obstimes), start=start_obs_id):
            irfs = self.irf_node.load_irfs(obstime)

            obs = Observation.create(
                obs_id=obs_id,
                pointing=self.pointing,
                livetime=float(obstime) * u.h,
                irfs=irfs,
                location=self.irf_node.location,
            )
            obs.meta["obstime"] = str(obstime)
            observations.append(obs)

        self.observations = Observations(observations)
        return self.observations

    def _resolve_observations(
        self,
        observations: Observations | None = None,
    ) -> Observations:
        if observations is not None:
            return observations
        if self.observations is not None:
            return self.observations
        raise ValueError("No Observations provided and no cached observations found")

    def create_spectrum_dataset_onoff(
        self,
        observation: Observation,
        *,
        dataset_name: str | None = None,
    ) -> SpectrumDatasetOnOff:
        energy_axis_reco = observation.bkg.axes["energy"]
        energy_axis_true = MapAxis.from_energy_bounds(
            0.3 * energy_axis_reco.edges[0],
            3.0 * energy_axis_reco.edges[-1],
            nbin=3 * len(energy_axis_reco.edges),
            name="energy_true",
        )

        geom = RegionGeom.create(region=self.on_region, axes=[energy_axis_reco])
        dataset_empty = SpectrumDataset.create(
            geom=geom,
            energy_axis_true=energy_axis_true,
            name=dataset_name or f"{self.source.name}_obs",
        )

        maker = SpectrumDatasetMaker(
            containment_correction=False,
            use_region_center=True,
            selection=["exposure", "edisp", "background"],
        )
        safe_mask_maker = SafeMaskMaker(methods=["bkg-peak"])

        dataset = maker.run(dataset_empty, observation)
        dataset = safe_mask_maker.run(dataset, observation)
        dataset.models = self.source.spectral_model.copy()

        return SpectrumDatasetOnOff.from_spectrum_dataset(
            dataset=dataset,
            acceptance=1,
            acceptance_off=self.n_off_regions,
        )

    def create_datasets_onoff(
        self,
        observations: Observations | None = None,
    ) -> dict[float, SpectrumDatasetOnOff]:
        observations = self._resolve_observations(observations)

        datasets_by_obstime: dict[float, SpectrumDatasetOnOff] = {}
        for observation in observations:
            obstime = float(observation.meta["obstime"])
            datasets_by_obstime[obstime] = self.create_spectrum_dataset_onoff(
                observation,
                dataset_name=f"{self.source.name}_{obstime:g}h",
            )

        self.datasets_onoff = dict(sorted(datasets_by_obstime.items()))
        return self.datasets_onoff

    def fake_data(
        self,
        dataset_on_off: SpectrumDatasetOnOff,
        *,
        n_fake_datasets: int | None = None,
    ) -> tuple[Datasets, QTable, np.ndarray]:
        n_fake = self.n_fake_datasets if n_fake_datasets is None else n_fake_datasets

        datasets = Datasets()
        for idx in range(n_fake):
            ds = dataset_on_off.copy(name=f"{dataset_on_off.name}_{idx}")
            # Dataset.copy() doesnt really work with its models - need to copy model by hand
            ds.models = self.source.spectral_model.copy()
            ds.fake(random_state=idx, npred_background=ds.npred_background())
            ds.meta_table["OBS_ID"] = [idx]
            datasets.append(ds)

        info_table = datasets.info_table()
        sigma = np.asarray(info_table["sqrt_ts"], dtype=float)

        return datasets, info_table, sigma

    def run_fake_studies(
        self,
        observations: Observations | None = None,
        *,
        n_fake_datasets: int | None = None,
    ) -> dict[float, dict[str, Any]]:
        datasets_onoff = self.create_datasets_onoff(observations)

        results: dict[float, dict[str, Any]] = {}
        for obstime, dataset_on_off in datasets_onoff.items():
            fake_datasets, info_table, sigma = self.fake_data(
                dataset_on_off,
                n_fake_datasets=n_fake_datasets,
            )

            results[obstime] = {
                "dataset_onoff": dataset_on_off,
                "fake_datasets": fake_datasets,
                "info_table": info_table,
                "sigma": sigma,
                "sigma_mean": float(np.mean(sigma)),
                "sigma_std": float(np.std(sigma)),
            }

        self.fake_results = dict(sorted(results.items()))
        return self.fake_results

    @staticmethod
    def _fit_obstime_from_significance(
        sigma_by_obstime: dict[float, np.ndarray],
        sigma_target: float,
    ) -> tuple[float, float]:
        obstimes = np.array(sorted(sigma_by_obstime), dtype=float)
        sigma_mean = np.array(
            [np.mean(sigma_by_obstime[t]) for t in obstimes],
            dtype=float,
        )
        sigma_std = np.array(
            [np.std(sigma_by_obstime[t]) for t in obstimes],
            dtype=float,
        )

        valid = np.isfinite(sigma_mean) & np.isfinite(sigma_std) & (sigma_std > 0)
        if valid.sum() < 2:
            return np.nan, np.nan

        # https://arxiv.org/pdf/2505.21632v1
        def model(t, a):
            return a * np.sqrt(t)

        popt, pcov = curve_fit(
            model,
            obstimes[valid],
            sigma_mean[valid],
            sigma=sigma_std[valid],
            absolute_sigma=True,
        )
        a = float(popt[0])
        da = float(np.sqrt(pcov[0, 0]))

        obstime_pred = (sigma_target / a) ** 2
        obstime_std = obstime_pred * 2.0 * (da / a)

        return obstime_pred, obstime_std

    def estimate_obstime_5sigma(
        self,
        observations: Observations | None = None,
        *,
        sigma_target: float = 5.0,
        n_fake_datasets: int | None = None,
        max_extrapolation_factor: float = 5.0,
    ) -> tuple[float, float, dict[float, dict[str, Any]]]:
        results = self.run_fake_studies(
            observations=observations,
            n_fake_datasets=n_fake_datasets,
        )

        sigma_by_obstime = {
            obstime: payload["sigma"] for obstime, payload in results.items()
        }

        obstime_5s, obstime_5s_std = self._fit_obstime_from_significance(
            sigma_by_obstime=sigma_by_obstime,
            sigma_target=sigma_target,
        )

        max_simulated_obstime = max(sigma_by_obstime)

        if np.isfinite(obstime_5s):
            if obstime_5s > max_extrapolation_factor * max_simulated_obstime:
                obstime_5s = np.nan
                obstime_5s_std = np.nan

        self.obstime_5s = obstime_5s
        self.obstime_5s_std = obstime_5s_std

        return obstime_5s, obstime_5s_std, results

    def plot_source_model_with_sensitivities(
        self,
        *,
        out_path: str | Path | None = None,
        energy_bounds: u.Quantity = [0.01, 100.0] * u.TeV,
    ):
        fig, ax = plt.subplots()
        e_lim = [5.0e-3, 5.0e2]

        for obstime in self.irf_node.common_obstimes:
            sens = self.irf_node.load_benchmark(obstime)
            energy_center = 0.5 * (sens["ENERG_LO"] + sens["ENERG_HI"])
            xerr = 0.5 * (sens["ENERG_HI"] - sens["ENERG_LO"])

            ax.errorbar(
                energy_center.flatten(),
                sens["ENERGY_FLUX_SENSITIVITY"].flatten(),
                xerr=xerr,
                ls="",
                label=f"CTAO-N - {obstime:g}h",
            )

        ax = self.source.spectral_model.spectral_model.plot(
            energy_bounds=energy_bounds,
            ax=ax,
            label=self.source.name,
            sed_type="e2dnde",
        )

        add_sensitivity_comparisons(ax, energy_limits=e_lim)

        ax.set_ylim(3.0e-14, 1.0e-9)
        ax.set_xlim(e_lim)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_title("Sensitivity")
        ax.set_xlabel(r"$E_{True}$ / TeV")
        ax.set_ylabel(
            r"$E^{2} \times$ Flux Sensitivity / $erg \cdot cm^{-2} \cdot s^{-1}$"
        )
        ax.grid(which="both", linestyle=":")
        ax.legend(loc="upper right", fontsize="small")
        fig.tight_layout()

        if out_path is not None:
            out_path = Path(out_path)
            out_path.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(out_path)

        return fig, ax

    def append_results_to_source(
        self,
        *,
        results: dict[float, dict[str, Any]] | None = None,
        sigma_target: float = 5.0,
    ) -> None:
        if results is None:
            if self.fake_results is None:
                raise ValueError("No fake-study results available")
            results = self.fake_results

        values: dict[str, Any] = {
            "matched_node_zen": self.irf_node.zen,
            "matched_node_az": self.irf_node.az,
            "matched_node_delta_b": self.irf_node.delta_b,
            "matched_node_sin_delta": self.irf_node.sin_delta,
            "matched_node_cos_theta": self.irf_node.cos_theta,
            "obstime_5s": self.obstime_5s,
            "obstime_5s_std": self.obstime_5s_std,
            "sigma_target": sigma_target,
        }

        for obstime, payload in sorted(results.items()):
            label = self._format_obstime(obstime)
            values[f"sigma_{label}h_mean"] = payload["sigma_mean"]
            values[f"sigma_{label}h_std"] = payload["sigma_std"]

        self.source.append_columns(values)


def main(source, output, irf_paths, benchmark_paths):
    source = Source(source)

    irfs = SubarrayIRFs(irf_paths=irf_paths, benchmark_paths=benchmark_paths)

    node = irfs.get_nearest_node(
        cos_theta_mean=source.row["cos_theta_mean"],
        sin_delta_mean=source.row["sin_delta_mean"],
    )

    analysis = SourceAnalysis(
        source=source,
        irf_node=node,
    )

    analysis.create_observations()
    _, _, results = analysis.estimate_obstime_5sigma()

    analysis.plot_source_model_with_sensitivities(
        out_path=Path(output).parent / f"{source.name}_sensitivity.pdf",
    )

    analysis.append_results_to_source(results=results)
    source.write(output)


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
