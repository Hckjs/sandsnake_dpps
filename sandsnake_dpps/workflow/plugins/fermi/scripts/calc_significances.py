import re
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path
from typing import Any, Iterable

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
from astropy.table import QTable
from astropy.time import Time
from gammapy.data import Observation
from gammapy.datasets import SpectrumDataset
from gammapy.irf import load_irf_dict_from_file
from gammapy.makers import SafeMaskMaker, SpectrumDatasetMaker
from gammapy.maps import MapAxis, RegionGeom
from gammapy.modeling.models import (
    EBLAbsorptionNormSpectralModel,
    ExpCutoffPowerLawNormSpectralModel,
    LogParabolaSpectralModel,
    PowerLawSpectralModel,
    SkyModel,
    SuperExpCutoffPowerLaw4FGLDR3SpectralModel,
)
from gammapy.stats import WStatCountsStatistic
from regions import PointSkyRegion

from common.plotting.colors import CTAO_COLORS
from core.scripts.mc.irf_plots import add_sensitivity_comparisons
from plugins.fermi.scripts.catalog_priors import RedshiftSource, SourceOrigin
from plugins.fermi.scripts.process_catalog import (
    CATALOG_NAMES,
    VisibilityConfig,
    get_B_direction,
)


site_params = VisibilityConfig()
REFERENCE_OBSTIME_H = 50.0


class AnalysisStatus(StrEnum):
    NOT_RUN = "not_run"
    SUCCESS = "success"
    EXTENDED_SOURCE = "extended_source"
    NOT_OBSERVABLE = "not_observable"
    UNKNOWN_ORIGIN = "unknown_origin"
    NO_USABLE_REDSHIFT = "no_usable_redshift"


class RedshiftScenarioLabel(StrEnum):
    LOW = "z_low"
    MED = "z_med"
    HIGH = "z_high"
    MEASURED = "z_measured"


REDSHIFT_PRIOR_SCENARIOS = (
    (RedshiftScenarioLabel.LOW, "z_q_low"),
    (RedshiftScenarioLabel.MED, "z_q_med"),
    (RedshiftScenarioLabel.HIGH, "z_q_high"),
)


@dataclass(frozen=True, slots=True)
class RedshiftScenario:
    label: str
    redshift: float
    source_column: str | None = None


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

    irf_path: Path
    benchmark_path: Path

    @property
    def key(self) -> tuple[int, int]:
        return (
            int(self.zen.to_value(u.deg)),
            int(self.az.to_value(u.deg)),
        )

    def load_benchmark(self) -> QTable:
        return QTable.read(self.benchmark_path, hdu="SENSITIVITY")

    def load_irfs(self):
        return load_irf_dict_from_file(self.irf_path)


class IRFCollection:
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
        b_declination: u.Quantity = site_params.prod_site_B_declination,
        b_inclination: u.Quantity = site_params.prod_site_B_inclination,
    ):
        self.site_name = site_name
        self.location = EarthLocation.of_site(site_name)
        self.prod_site_B_declination = b_declination
        self.prod_site_B_inclination = b_inclination
        self.frame = AltAz(
            location=self.location, obstime=Time("2027-01-01T00:00:00", scale="utc")
        )
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

    def _paths_by_node(
        self,
        paths: Iterable[str | Path],
    ) -> dict[tuple[int, int], Path]:
        grouped: dict[tuple[int, int], Path] = {}

        for path in paths:
            zen, az, obstime, parsed_path = self.parse_input_path(path)
            key = (int(zen.to_value(u.deg)), int(az.to_value(u.deg)))
            if obstime != REFERENCE_OBSTIME_H:
                raise ValueError(
                    f"Expected only {REFERENCE_OBSTIME_H:g} h inputs, got {parsed_path}"
                )
            if key in grouped:
                raise ValueError(
                    f"Duplicate input for node {key}: {grouped[key]} and {parsed_path}"
                )
            grouped[key] = parsed_path

        return dict(grouped)

    def create_node(
        self,
        irf_path: Path,
        benchmark_path: Path,
    ) -> IRFNode:
        zen_ref, az_ref, _, _ = self.parse_input_path(irf_path)
        for path in (irf_path, benchmark_path):
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
            irf_path=irf_path,
            benchmark_path=benchmark_path,
        )

    def build_nodes(
        self,
        irf_paths: Iterable[str | Path],
        benchmark_paths: Iterable[str | Path],
    ) -> dict[tuple[int, int], IRFNode]:
        grouped_irfs = self._paths_by_node(irf_paths)
        grouped_benchmarks = self._paths_by_node(benchmark_paths)
        if grouped_irfs.keys() != grouped_benchmarks.keys():
            raise ValueError("50 h IRF and benchmark nodes do not match")

        return {
            key: self.create_node(grouped_irfs[key], grouped_benchmarks[key])
            for key in sorted(grouped_irfs)
        }

    def get_nearest_node_key(
        self,
        cos_theta_mean: float,
        sin_delta_mean: float,
        weight_cos_theta: float = 1.0,
        weight_sin_delta: float = 1.0,
    ) -> tuple[int, int]:
        if not self.nodes:
            raise ValueError("No nodes available in IRFCollection")

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

    This class owns the catalog-derived source properties, redshift-scenario
    resolution, and the corresponding spectral models. It does not own IRF,
    dataset or significance-estimation state.
    """

    def __init__(
        self,
        path: str | Path,
        *,
        ebl_model: str = "saldana-lopez21",
    ):
        self.path = Path(path)
        self.table = QTable.read(self.path, format="ascii.ecsv")

        if len(self.table) != 1:
            raise ValueError(
                f"Expected exactly one row in source table {self.path}, "
                f"got {len(self.table)}"
            )

        self.ebl_model = ebl_model

        if "catalog" not in self.table.colnames:
            raise ValueError(
                f"Source table {self.path} is missing required column 'catalog'"
            )
        self.catalog = str(self.row["catalog"])
        valid_catalogs = {
            CATALOG_NAMES["4FGL"],
            CATALOG_NAMES["3FHL"],
            CATALOG_NAMES["4FHL"],
        }
        if self.catalog not in valid_catalogs:
            raise ValueError(
                f"Unknown catalog {self.catalog!r} in {self.path}; "
                f"expected one of {sorted(valid_catalogs)}"
            )

        self.position = SkyCoord(
            self.row["RAJ2000"],
            self.row["DEJ2000"],
            unit="deg",
            frame="icrs",
        )
        self.origin = self.row["z_class"]
        self.z_source = self.row["z_source"]
        ext_value = self.row["is_extended_source"]
        self.is_extended_source = (
            False if np.ma.is_masked(ext_value) else bool(ext_value)
        )

        if self.is_extended_source:
            self.redshift_scenarios = None
            self.spectral_models = None
        else:
            self.redshift_scenarios = self._resolve_redshift_scenarios()
            self.spectral_models = self._create_spectral_model()

    @property
    def row(self):
        return self.table[0]

    @property
    def name(self) -> str:
        return str(self.row["Source_Name"])

    @property
    def redshift_by_label(self) -> dict[str, float] | None:
        if self.redshift_scenarios is None:
            return None

        return {
            scenario.label: scenario.redshift for scenario in self.redshift_scenarios
        }

    @property
    def spectral_model_by_label(self) -> dict[str, SkyModel] | None:
        if self.spectral_models is None:
            return None

        return {model.name: model for model in self.spectral_models}

    @property
    def has_redshift_scenarios(self) -> bool:
        return bool(self.redshift_scenarios)

    @property
    def has_prior_redshift_scenarios(self) -> bool:
        return self.z_source in [RedshiftSource.PRIOR_BLL, RedshiftSource.PRIOR_FSRQ]

    def _resolve_redshift_scenarios(self) -> list[RedshiftScenario] | None:
        if self.has_prior_redshift_scenarios:
            return [
                RedshiftScenario(
                    label=label,
                    redshift=self.row[column],
                    source_column=column,
                )
                for label, column in REDSHIFT_PRIOR_SCENARIOS
            ]

        if self.z_source == RedshiftSource.MEASURED:
            return [
                RedshiftScenario(
                    label=RedshiftScenarioLabel.MEASURED,
                    redshift=self.row[RedshiftScenarioLabel.MEASURED],
                    source_column=RedshiftScenarioLabel.MEASURED,
                )
            ]

        return None

    def _create_base_spectral_model(self):
        if self.catalog == CATALOG_NAMES["3FHL"]:
            spec_type = str(self.row["SpectrumType"])
            if spec_type == "PowerLaw":
                return PowerLawSpectralModel(
                    amplitude=self.row["Flux_Density"] / u.ph,
                    reference=self.row["Pivot_Energy"],
                    index=self.row["PowerLaw_Index"],
                )

            if spec_type == "LogParabola":
                return LogParabolaSpectralModel(
                    amplitude=self.row["Flux_Density"] / u.ph,
                    reference=self.row["Pivot_Energy"],
                    alpha=self.row["Spectral_Index"],
                    beta=self.row["beta"],
                )

            raise ValueError(
                f"Spectral model {spec_type!r} not implemented for 3FHL source"
            )

        if self.catalog == CATALOG_NAMES["4FGL"]:
            spec_type = str(self.row["SpectrumType"])
            if spec_type == "PowerLaw":
                return PowerLawSpectralModel(
                    amplitude=self.row["PL_Flux_Density"] / u.ph,
                    reference=self.row["Pivot_Energy"],
                    index=self.row["PL_Index"],
                )

            if spec_type == "LogParabola" or (
                spec_type == "PLSuperExpCutoff"
                and self.origin == SourceOrigin.EXTRAGALACTIC
            ):
                return LogParabolaSpectralModel(
                    amplitude=self.row["LP_Flux_Density"] / u.ph,
                    reference=self.row["Pivot_Energy"],
                    alpha=self.row["LP_Index"],
                    beta=self.row["LP_beta"],
                )

            if spec_type == "PLSuperExpCutoff":
                return SuperExpCutoffPowerLaw4FGLDR3SpectralModel(
                    amplitude=self.row["PLEC_Flux_Density"] / u.ph,
                    reference=self.row["Pivot_Energy"],
                    index_1=self.row["PLEC_IndexS"],
                    index_2=self.row["PLEC_Exp_Index"],
                    expfactor=self.row["PLEC_ExpfactorS"],
                )

            raise ValueError(f"Spectral model {spec_type!r} not implemented")

        if self.catalog == CATALOG_NAMES["4FHL"]:
            reference = 100 * u.GeV
            energy_min = 50 * u.GeV
            energy_max = 2 * u.TeV
            flux50 = u.Quantity(self.row["Flux50"], copy=False)
            if flux50.unit == u.dimensionless_unscaled:
                flux50 = flux50 * u.Unit("cm-2 s-1")
            flux50 = flux50.to("cm-2 s-1")
            if not np.isfinite(flux50.value) or flux50 <= 0 * flux50.unit:
                raise ValueError(f"Invalid 4FHL Flux50 for {self.name}: {flux50}")

            model = PowerLawSpectralModel(
                amplitude=1 * u.Unit("cm-2 s-1 TeV-1"),
                reference=reference,
                index=float(self.row["PL_Index"]),
            )
            scale = (flux50 / model.integral(energy_min, energy_max)).to_value("")
            model.amplitude.value *= scale
            integrated_flux = model.integral(energy_min, energy_max)
            if not u.isclose(integrated_flux, flux50, rtol=1e-10):
                raise ValueError(
                    f"4FHL power-law normalization failed for {self.name}: "
                    f"integral={integrated_flux}, Flux50={flux50}"
                )
            return model

        raise ValueError(
            f"Unknown catalog {self.catalog!r} for spectral model creation"
        )

    def _create_cutoff_model(self, redshift: float):
        cutoff_energy = 10 * u.TeV / (1.0 + redshift)

        return ExpCutoffPowerLawNormSpectralModel(
            norm=1.0,
            index=0.0,
            lambda_=1.0 / cutoff_energy,
            alpha=1.0,
            reference=cutoff_energy / 10.0,
        )

    def _create_spectral_model(self) -> list[SkyModel] | None:
        base_model = self._create_base_spectral_model()
        if self.origin == SourceOrigin.GALACTIC:
            return [SkyModel(spectral_model=base_model, name=SourceOrigin.GALACTIC)]

        if self.redshift_scenarios is None:
            return None

        spec_model = []
        for scenario in self.redshift_scenarios:
            ebl_model = EBLAbsorptionNormSpectralModel.read_builtin(
                self.ebl_model,
                redshift=scenario.redshift,
            )
            cutoff_model = self._create_cutoff_model(scenario.redshift)
            spec_model.append(
                SkyModel(
                    spectral_model=base_model * ebl_model * cutoff_model,
                    name=scenario.label,
                )
            )

        return spec_model


class SourceAnalysis:
    """Analyze one source using the fixed, 50 h response of one IRF node."""

    def __init__(
        self,
        source: Source,
        irf_node: IRFNode,
        output_table: QTable,
        *,
        offset: u.Quantity | None = None,  # Needs to match IRFs
        n_off_regions: int = 5,
    ):
        self.source = source
        self.irf_node = irf_node
        self.output_table = output_table
        self.offset = offset
        self.n_off_regions = n_off_regions

        if offset is not None and offset.value > 0.0:
            self.pointing = self.source.position.directional_offset_by(
                position_angle=90 * u.deg,
                separation=self.offset,
            )
        else:
            self.pointing = self.source.position

        self.on_region = PointSkyRegion(self.source.position)

    @property
    def alpha(self) -> float:
        return 1.0 / self.n_off_regions

    def run(self) -> None:
        if self.source.spectral_models is None:
            raise ValueError("Source has no spectral model(s) to run analysis with")

        self.output_table["matched_node_zen"] = [self.irf_node.zen]
        self.output_table["matched_node_az"] = [self.irf_node.az]
        self.output_table["matched_node_delta_b"] = [self.irf_node.delta_b]
        self.output_table["matched_node_sin_delta"] = [self.irf_node.sin_delta]
        self.output_table["matched_node_cos_theta"] = [self.irf_node.cos_theta]

        observation = self._create_observation()
        for spectral_model in self.source.spectral_models:
            dataset = self.create_spectrum_dataset(observation, spectral_model)
            sigma_asimov = self.compute_asimov_significance(dataset)
            obstime = self.estimate_obstime_from_reference(
                sigma_asimov=sigma_asimov,
                reference_obstime=REFERENCE_OBSTIME_H,
                sigma_target=float(self.output_table[0]["sigma_target"]),
            )
            self.append_results_to_output_table(
                sigma_asimov, obstime, spectral_model.name
            )

        self.output_table["status"] = [AnalysisStatus.SUCCESS]

    def _create_observation(self) -> Observation:
        return Observation.create(
            obs_id=1,
            pointing=self.pointing,
            livetime=REFERENCE_OBSTIME_H * u.h,
            irfs=self.irf_node.load_irfs(),
            location=self.irf_node.location,
        )

    def _background_counts_from_rad_max(
        self,
        observation: Observation,
        energy_axis_reco: MapAxis,
    ) -> np.ndarray:
        if observation.rad_max is None:
            raise ValueError(
                "Point-like RAD_MAX background requested, but observation.rad_max "
                "is missing."
            )
        if observation.bkg is None:
            raise ValueError(
                "Point-like RAD_MAX background requested, but observation.bkg "
                "is missing."
            )

        source_offset = 0.0 * u.deg if self.offset is None else self.offset
        energy = energy_axis_reco.center
        energy_width = np.diff(energy_axis_reco.edges)
        rad_max = observation.rad_max.evaluate(energy=energy, offset=source_offset)
        theta = rad_max.to(u.rad)
        solid_angle = 2.0 * np.pi * (1.0 - np.cos(theta.value)) * u.sr
        bkg_rate = observation.bkg.evaluate(energy=energy, offset=source_offset)
        livetime = observation.observation_live_time_duration
        background_counts = (
            bkg_rate * energy_width * livetime * solid_angle
        ).to_value("")
        background_counts = np.asarray(background_counts, dtype=float).reshape(-1)
        background_counts[~np.isfinite(background_counts)] = 0.0
        background_counts[background_counts < 0.0] = 0.0

        if background_counts.size != energy_axis_reco.nbin:
            raise ValueError(
                "RAD_MAX background shape mismatch: "
                f"got {background_counts.size} bins, expected {energy_axis_reco.nbin}"
            )
        if not np.any(background_counts > 0.0):
            raise ValueError(
                "RAD_MAX background is zero in all reconstructed-energy bins. "
                f"source_offset={source_offset:.3f}"
            )
        return background_counts

    def _set_rad_max_background(
        self,
        dataset: SpectrumDataset,
        observation: Observation,
        energy_axis_reco: MapAxis,
    ) -> None:
        background_counts = self._background_counts_from_rad_max(
            observation, energy_axis_reco
        )
        background_data = np.zeros(dataset.counts.data.shape, dtype=float)
        background_data[...] = background_counts.reshape(
            (energy_axis_reco.nbin,) + (1,) * (background_data.ndim - 1)
        )
        dataset.background = dataset.counts.copy(data=background_data)

    def create_spectrum_dataset(
        self,
        observation: Observation,
        spectral_model: SkyModel,
    ) -> SpectrumDataset:
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
            name=f"{self.source.name}_{spectral_model.name}_{REFERENCE_OBSTIME_H:g}h",
        )
        dataset = SpectrumDatasetMaker(
            containment_correction=False, selection=["exposure", "edisp"]
        ).run(dataset_empty, observation)
        self._set_rad_max_background(dataset, observation, energy_axis_reco)
        dataset = SafeMaskMaker(methods=["aeff-default"]).run(dataset, observation)
        dataset.models = spectral_model.copy()
        return dataset

    def compute_asimov_significance(self, dataset: SpectrumDataset) -> float:
        """Return the expected Li & Ma ON/OFF excess significance."""
        signal = np.asarray(dataset.npred_signal().data, dtype=float)
        background = np.asarray(dataset.npred_background().data, dtype=float)
        if dataset.mask_safe is None:
            safe = np.ones(signal.shape, dtype=bool)
        else:
            safe = np.broadcast_to(
                np.asarray(dataset.mask_safe.data, dtype=bool), signal.shape
            )

        s = float(np.sum(signal[safe]))
        b = float(np.sum(background[safe]))
        stat = WStatCountsStatistic(
            n_on=s + b,
            n_off=b / self.alpha,
            alpha=self.alpha,
        )
        return float(np.asarray(stat.sqrt_ts))

    @staticmethod
    def estimate_obstime_from_reference(
        sigma_asimov: float,
        reference_obstime: float,
        sigma_target: float,
    ) -> float:
        """Scale a fixed-response Asimov significance to the target time.

        Signal and background expectations scale linearly with time while alpha
        stays constant.  The 50 h IRF, including its 50 h optimized cuts, is held
        fixed at every scaled time, so S_A(T) is proportional to sqrt(T).  The
        result is therefore not based on cuts re-optimized for the predicted time.
        """
        if not np.isfinite(sigma_asimov) or sigma_asimov <= 0:
            return np.nan
        return float(reference_obstime * (sigma_target / sigma_asimov) ** 2)

    def append_results_to_output_table(
        self, sigma_asimov: float, obstime_target: float, model_label: str
    ) -> None:
        suffix = f"_{model_label}" if self.source.has_prior_redshift_scenarios else ""
        self.output_table[f"sigma_asimov_50h{suffix}"] = sigma_asimov
        self.output_table[f"obstime_5s{suffix}"] = obstime_target

    def write(self, outpath: Path, *, overwrite: bool = True) -> None:
        self.output_table.write(outpath, format="ascii.ecsv", overwrite=overwrite)

    @staticmethod
    def _sed_flux(model: SkyModel, energy: u.Quantity) -> u.Quantity:
        return (energy**2 * model.spectral_model(energy)).to("erg cm-2 s-1")

    @staticmethod
    def _row_quantity(row, column: str, default_unit: str | u.Unit) -> u.Quantity:
        value = row[column]
        quantity = u.Quantity(value, copy=False)

        values = np.ma.filled(quantity.value, np.nan)
        unit = quantity.unit
        default_unit = u.Unit(default_unit)

        if unit == u.dimensionless_unscaled:
            unit = default_unit

        return u.Quantity(values, unit).to(default_unit)

    def _fermi_flux_point_energy_edges(self, n_bins: int) -> u.Quantity | None:
        """
        Return catalog SED-bin edges.

        4FGL-DR4:
            50 MeV - 1 TeV, 8 bins.

        3FHL:
            10 GeV - 2 TeV, 5 bins.

        4FHL:
            50 GeV - 2 TeV, 3 bins.
        """
        catalog = self.source.catalog

        if catalog == CATALOG_NAMES["3FHL"]:
            edges = [10, 20, 50, 150, 500, 2000] * u.GeV

        elif catalog == CATALOG_NAMES["4FHL"]:
            edges = [50, 171, 585, 2000] * u.GeV

        elif catalog == CATALOG_NAMES["4FGL"]:
            edges = [
                50,
                100,
                300,
                1_000,
                3_000,
                10_000,
                30_000,
                100_000,
                1_000_000,
            ] * u.MeV

        else:
            raise ValueError(f"Unknown catalog {catalog!r} for Fermi flux-point bins")

        if len(edges) != n_bins + 1:
            log.warning(
                "%s: expected %d flux-point edges for %d bins, got %d edges",
                self.source.name,
                n_bins + 1,
                n_bins,
                len(edges),
            )
            return None

        return edges.to(u.TeV)

    def _fermi_flux_points_from_row(self) -> dict[str, Any] | None:
        """
        Extract Fermi SED points from the catalog row.

        Returns e2dnde-like points in erg cm-2 s-1, including asymmetric
        uncertainties and upper-limit values.
        """
        row = self.source.row

        if self.source.catalog == CATALOG_NAMES["4FHL"]:
            return self._fermi_flux_points_from_4fhl_row()

        required = {"Flux_Band", "Unc_Flux_Band", "Sqrt_TS_Band"}
        if not required.issubset(row.colnames):
            return None

        if "nuFnu_Band" in row.colnames:
            nufnu_column = "nuFnu_Band"
        elif "nuFnu" in row.colnames:
            nufnu_column = "nuFnu"
        else:
            return None

        flux = self._row_quantity(row, "Flux_Band", "ph cm-2 s-1")
        flux_err = self._row_quantity(row, "Unc_Flux_Band", "ph cm-2 s-1")
        e2dnde = self._row_quantity(row, nufnu_column, "erg cm-2 s-1")
        sqrt_ts = np.asarray(np.ma.filled(row["Sqrt_TS_Band"], np.nan), dtype=float)

        n_bins = len(e2dnde)
        edges = self._fermi_flux_point_energy_edges(n_bins)
        if edges is None:
            return None

        e_min = edges[:-1]
        e_max = edges[1:]
        e_ref = np.sqrt(e_min * e_max)
        xerr = u.Quantity(
            [
                (e_ref - e_min).to_value(u.TeV),
                (e_max - e_ref).to_value(u.TeV),
            ],
            u.TeV,
        )

        flux_value = flux.to_value("ph cm-2 s-1")
        flux_err_value = flux_err.to_value("ph cm-2 s-1")
        e2dnde_value = e2dnde.to_value("erg cm-2 s-1")

        with np.errstate(divide="ignore", invalid="ignore"):
            e2dnde_errn_value = np.abs(e2dnde_value * flux_err_value[:, 0] / flux_value)
            e2dnde_errp_value = e2dnde_value * flux_err_value[:, 1] / flux_value

        e2dnde_errn = e2dnde_errn_value * u.Unit("erg cm-2 s-1")
        e2dnde_errp = e2dnde_errp_value * u.Unit("erg cm-2 s-1")

        # Fermi/Gammapy convention:
        # lower error NaN -> upper limit.
        # Sqrt_TS_Band < 1 is also treated as an upper limit.
        is_ul = ~np.isfinite(e2dnde_errn_value) | (sqrt_ts < 1.0)

        e2dnde_ul = e2dnde + 2.0 * e2dnde_errp
        invalid_ul = ~np.isfinite(e2dnde_ul.to_value("erg cm-2 s-1"))
        if np.any(invalid_ul):
            e2dnde_ul[invalid_ul] = e2dnde[invalid_ul]

        catalog_label = {
            CATALOG_NAMES["4FGL"]: "4FGL-DR4",
            CATALOG_NAMES["3FHL"]: "3FHL",
        }[self.source.catalog]

        return {
            "catalog_label": catalog_label,
            "e_ref": e_ref,
            "xerr": xerr,
            "e2dnde": e2dnde,
            "e2dnde_errn": e2dnde_errn,
            "e2dnde_errp": e2dnde_errp,
            "e2dnde_ul": e2dnde_ul,
            "sqrt_ts": sqrt_ts,
            "is_ul": is_ul,
        }

    def _fermi_flux_points_from_4fhl_row(self) -> dict[str, Any] | None:
        """Convert the three 4FHL integral-flux bands to differential SED points."""
        row = self.source.row
        flux_columns = (
            "Flux50_171GeV",
            "Flux171_585GeV",
            "Flux585_2000GeV",
        )
        error_columns = (
            "Unc_Flux50_171GeV",
            "Unc_Flux171_585GeV",
            "Unc_Flux585_2000GeV",
        )
        sqrt_ts_columns = (
            "Sqrt_TS50_171GeV",
            "Sqrt_TS171_585GeV",
            "Sqrt_TS585_2000GeV",
        )
        required = set(flux_columns + error_columns + sqrt_ts_columns)
        if not required.issubset(row.colnames):
            return None

        flux = u.Quantity(
            [u.Quantity(row[column], copy=False).value for column in flux_columns],
            u.Unit("cm-2 s-1"),
        )
        flux_err = u.Quantity(
            [u.Quantity(row[column], copy=False).value for column in error_columns],
            u.Unit("cm-2 s-1"),
        )
        sqrt_ts = np.asarray([row[column] for column in sqrt_ts_columns], dtype=float)
        edges = self._fermi_flux_point_energy_edges(len(flux))
        e_min, e_max = edges[:-1], edges[1:]
        e_ref = np.sqrt(e_min * e_max)

        unit_model = PowerLawSpectralModel(
            amplitude=1 * u.Unit("cm-2 s-1 TeV-1"),
            reference=100 * u.GeV,
            index=float(row["PL_Index"]),
        )
        conversion = (
            e_ref**2 * unit_model(e_ref) / unit_model.integral(e_min, e_max)
        ).to(u.erg)
        e2dnde = (flux * conversion).to("erg cm-2 s-1")
        e2dnde_error = (flux_err * conversion).to("erg cm-2 s-1")
        e2dnde_errn = e2dnde_error.copy()
        e2dnde_errp = e2dnde_error.copy()
        is_ul = (
            ~np.isfinite(e2dnde.value)
            | ~np.isfinite(e2dnde_errn.value)
            | (sqrt_ts < 1.0)
        )
        e2dnde_ul = e2dnde + 2.0 * e2dnde_errp
        xerr = u.Quantity(
            [
                (e_ref - e_min).to_value(u.TeV),
                (e_max - e_ref).to_value(u.TeV),
            ],
            u.TeV,
        )
        return {
            "catalog_label": "4FHL",
            "e_ref": e_ref,
            "xerr": xerr,
            "e2dnde": e2dnde,
            "e2dnde_errn": e2dnde_errn,
            "e2dnde_errp": e2dnde_errp,
            "e2dnde_ul": e2dnde_ul,
            "sqrt_ts": sqrt_ts,
            "is_ul": is_ul,
        }

    def _plot_fermi_flux_points(self, ax) -> None:
        flux_points = self._fermi_flux_points_from_row()
        if flux_points is None:
            log.debug("%s: no Fermi flux points found in source row", self.source.name)
            return

        e_ref = flux_points["e_ref"]
        xerr = flux_points["xerr"]
        y = flux_points["e2dnde"]
        yerrn = flux_points["e2dnde_errn"]
        yerrp = flux_points["e2dnde_errp"]
        y_ul = flux_points["e2dnde_ul"]
        is_ul = flux_points["is_ul"]
        catalog_label = flux_points["catalog_label"]

        y_value = y.to_value("erg cm-2 s-1")
        yerrn_value = yerrn.to_value("erg cm-2 s-1")
        yerrp_value = yerrp.to_value("erg cm-2 s-1")
        y_ul_value = y_ul.to_value("erg cm-2 s-1")

        is_point = (
            ~is_ul
            & np.isfinite(y_value)
            & np.isfinite(yerrn_value)
            & np.isfinite(yerrp_value)
            & (y_value > 0.0)
        )
        is_upper_limit = is_ul & np.isfinite(y_ul_value) & (y_ul_value > 0.0)

        if np.any(is_point):
            ax.errorbar(
                e_ref[is_point].to_value(u.TeV),
                y[is_point].to_value("erg cm-2 s-1"),
                xerr=[
                    xerr[0][is_point].to_value(u.TeV),
                    xerr[1][is_point].to_value(u.TeV),
                ],
                yerr=[
                    yerrn[is_point].to_value("erg cm-2 s-1"),
                    yerrp[is_point].to_value("erg cm-2 s-1"),
                ],
                fmt="o",
                ls="",
                markersize=4,
                capsize=2,
                label=f"{catalog_label} flux points",
                zorder=5,
                color=CTAO_COLORS["cherenkov_cyan"],
            )

        if np.any(is_upper_limit):
            # Matplotlib needs a finite yerr to draw the upper-limit arrow.
            # The point itself is placed at the upper-limit value.
            ul_yerr = 0.35 * y_ul[is_upper_limit].to_value("erg cm-2 s-1")

            ax.errorbar(
                e_ref[is_upper_limit].to_value(u.TeV),
                y_ul[is_upper_limit].to_value("erg cm-2 s-1"),
                xerr=[
                    xerr[0][is_upper_limit].to_value(u.TeV),
                    xerr[1][is_upper_limit].to_value(u.TeV),
                ],
                yerr=ul_yerr,
                uplims=True,
                fmt="v",
                ls="",
                markersize=4,
                capsize=2,
                label=f"{catalog_label} upper limits",
                zorder=5,
                color=CTAO_COLORS["cherenkov_cyan"],
                alpha=0.5,
            )

    def plot_source_model_with_sensitivities(
        self,
        *,
        out_path: str | Path | None = None,
        energy_bounds: u.Quantity | None = None,
    ):
        if self.source.spectral_models is None:
            raise ValueError("Source has no spectral model(s) to plot")

        if self.source.catalog == CATALOG_NAMES["4FGL"]:
            e_lim = [5.0e-5, 1.0e3]
            if energy_bounds is None:
                energy_bounds = [5.0e-5, 100.0] * u.TeV
        elif self.source.catalog in {
            CATALOG_NAMES["3FHL"],
            CATALOG_NAMES["4FHL"],
        }:
            e_lim = [5.0e-3, 1.0e3]
            if energy_bounds is None:
                energy_bounds = [5.0e-3, 100.0] * u.TeV

        fig, ax = plt.subplots()

        sens = self.irf_node.load_benchmark()
        energy_center = 0.5 * (sens["ENERG_LO"] + sens["ENERG_HI"])
        xerr = 0.5 * (sens["ENERG_HI"] - sens["ENERG_LO"])
        ax.errorbar(
            energy_center.flatten(),
            sens["ENERGY_FLUX_SENSITIVITY"].flatten(),
            xerr=xerr.flatten(),
            ls="",
            color=CTAO_COLORS["interstellar_indigo"],
            label=f"CTAO-N - {REFERENCE_OBSTIME_H:g}h",
        )
        self._plot_fermi_flux_points(ax)

        if self.source.has_prior_redshift_scenarios:
            energy = (
                np.geomspace(
                    energy_bounds[0].to_value(u.TeV),
                    energy_bounds[1].to_value(u.TeV),
                    256,
                )
                * u.TeV
            )
            y_low = self._sed_flux(
                self.source.spectral_model_by_label[RedshiftScenarioLabel.LOW],
                energy,
            )
            y_high = self._sed_flux(
                self.source.spectral_model_by_label[RedshiftScenarioLabel.HIGH],
                energy,
            )

            ax.fill_between(
                energy.to_value(u.TeV),
                np.minimum(y_low.value, y_high.value),
                np.maximum(y_low.value, y_high.value),
                color=CTAO_COLORS["cherenkov_cyan"],
                alpha=0.15,
            )

            for label, linestyle in [
                (RedshiftScenarioLabel.LOW, ":"),
                (RedshiftScenarioLabel.MED, "-"),
                (RedshiftScenarioLabel.HIGH, "--"),
            ]:
                z = self.source.redshift_by_label[label]
                self.source.spectral_model_by_label[label].spectral_model.plot(
                    energy_bounds=energy_bounds,
                    ax=ax,
                    label=f"{label}={z:.3g}",
                    sed_type="e2dnde",
                    color=CTAO_COLORS["cherenkov_cyan"],
                    linestyle=linestyle,
                )
        else:
            self.source.spectral_models[0].spectral_model.plot(
                energy_bounds=energy_bounds,
                ax=ax,
                label=self.source.name,
                sed_type="e2dnde",
                color=CTAO_COLORS["cherenkov_cyan"],
            )

        add_sensitivity_comparisons(ax, energy_limits=e_lim, add_prod5=False)

        ax.set_ylim(3.0e-14, 1.0e-9)
        ax.set_xlim(e_lim)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_title(self.source.name)
        ax.set_xlabel(r"$E_{True}$ / TeV")
        ax.set_ylabel(
            r"$E^{2} \times$ Flux Sensitivity / $erg \cdot cm^{-2} \cdot s^{-1}$"
        )
        ax.grid(which="both", linestyle=":")
        ax.legend(loc="upper right", fontsize="x-small")

        if out_path is not None:
            fig.savefig(out_path)

        return fig, ax


def normalize_label(value) -> str:
    if np.ma.is_masked(value):
        return ""

    if isinstance(value, bytes):
        value = value.decode()

    text = str(value).strip().lower()
    if text in {"", "--", "nan", "none", "masked"}:
        return ""

    return text


def source_validity_check(source: Source, output_table: QTable) -> bool:
    if source.is_extended_source:
        output_table["status"] = [AnalysisStatus.EXTENDED_SOURCE]
        return False

    if not (
        np.isfinite(source.row["cos_theta_mean"])
        and np.isfinite(source.row["sin_delta_mean"])
    ):
        output_table["status"] = [AnalysisStatus.NOT_OBSERVABLE]
        return False

    if source.origin == SourceOrigin.UNKNOWN:
        output_table["status"] = [AnalysisStatus.UNKNOWN_ORIGIN]
        return False

    if (
        not source.has_redshift_scenarios
        and source.origin == SourceOrigin.EXTRAGALACTIC
    ):
        output_table["status"] = [AnalysisStatus.NO_USABLE_REDSHIFT]
        return False

    return True


def create_output_table(
    source_table: QTable,
    *,
    sigma_target: float = 5.0,
) -> QTable:
    output_table = source_table.copy()
    output_table["status"] = [AnalysisStatus.NOT_RUN]
    output_table["sigma_target"] = [sigma_target]
    output_table["matched_node_zen"] = [np.nan * u.deg]
    output_table["matched_node_az"] = [np.nan * u.deg]
    output_table["matched_node_delta_b"] = [np.nan * u.deg]
    output_table["matched_node_sin_delta"] = [np.nan]
    output_table["matched_node_cos_theta"] = [np.nan]

    for label, _ in REDSHIFT_PRIOR_SCENARIOS:
        output_table[f"sigma_asimov_50h_{label}"] = [np.nan]
        output_table[f"obstime_5s_{label}"] = [np.nan]
    output_table["sigma_asimov_50h"] = [np.nan]
    output_table["obstime_5s"] = [np.nan]
    return output_table


def main(
    source_path: str | Path,
    output: str | Path,
    irf_paths: Iterable[str | Path],
    benchmark_paths: Iterable[str | Path],
    *,
    sigma_target: float = 5.0,
) -> None:
    outpath = Path(output)
    outpath.parent.mkdir(parents=True, exist_ok=True)

    source = Source(source_path)
    irf_collection = IRFCollection(irf_paths=irf_paths, benchmark_paths=benchmark_paths)
    output_table = create_output_table(source.table, sigma_target=sigma_target)

    if not source_validity_check(source, output_table):
        output_table.write(outpath, format="ascii.ecsv", overwrite=True)
        return

    nearest_node = irf_collection.get_nearest_node(
        cos_theta_mean=source.row["cos_theta_mean"],
        sin_delta_mean=source.row["sin_delta_mean"],
    )

    analysis = SourceAnalysis(
        source=source,
        irf_node=nearest_node,
        output_table=output_table,
    )
    analysis.run()

    analysis.plot_source_model_with_sensitivities(
        out_path=outpath.parent / f"{source.name}_sensitivity.pdf",
    )

    analysis.write(outpath)


def parse_args():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--source", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--irfs", nargs="+", required=True)
    parser.add_argument("--benchmarks", nargs="+", required=True)
    parser.add_argument("--sigma-target", type=float, default=5.0)
    return parser.parse_args()


def main_from_snakemake(snakemake):
    main(
        snakemake.input.source,
        snakemake.output[0],
        snakemake.input.irfs,
        snakemake.input.benchmarks,
    )


def main_from_args(args):
    main(
        args.source,
        args.output,
        args.irfs,
        args.benchmarks,
        sigma_target=args.sigma_target,
    )


if "snakemake" in globals():
    main_from_snakemake(snakemake)  # noqa: F821
elif __name__ == "__main__":
    args = parse_args()
    main_from_args(args)
