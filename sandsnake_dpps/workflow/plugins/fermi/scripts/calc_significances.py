from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable
import re
import logging

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
from astropy.table import QTable
from astropy.time import Time
from gammapy.data import Observation, Observations
from gammapy.datasets import Datasets, SpectrumDataset, SpectrumDatasetOnOff
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
from regions import PointSkyRegion
from scipy.optimize import curve_fit

from core.scripts.mc.irf_plots import add_sensitivity_comparisons
from common.plotting.colors import CTAO_COLORS
from plugins.fermi.scripts.process_catalog import VisibilityConfig, get_B_direction
from plugins.fermi.scripts.catalog_priors import RedshiftSource, SourceOrigin
from enum import StrEnum


site_params = VisibilityConfig()
log = logging.getLogger(__name__)


class AnalysisStatus(StrEnum):
    NOT_RUN = "not_run"
    SUCCESS = "success"
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

    obstime_paths_irfs: dict[float, Path]
    obstime_paths_benchmarks: dict[float, Path]

    @property
    def key(self) -> tuple[int, int]:
        return (
            int(self.zen.to_value(u.deg)),
            int(self.az.to_value(u.deg)),
        )

    @property
    def obstimes(self) -> list[float]:
        return sorted(self.obstime_paths_irfs)

    def get_irf_path(self, obstime: float) -> Path:
        return self.obstime_paths_irfs[float(obstime)]

    def get_benchmark_path(self, obstime: float) -> Path:
        return self.obstime_paths_benchmarks[float(obstime)]

    def load_benchmark(self, obstime: float) -> QTable:
        return QTable.read(self.get_benchmark_path(obstime), hdu="SENSITIVITY")

    def load_irfs(self, obstime: float):
        return load_irf_dict_from_file(self.get_irf_path(obstime))


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

    @property
    def obstimes(self) -> list[float]:
        obstimes: set[float] = set()
        for node in self.nodes.values():
            obstimes.update(node.obstimes)
        return sorted(obstimes)

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

        for path in {**irf_paths, **benchmark_paths}.values():
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
    dataset, fake-data, or significance-estimation state.
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

        self.position = SkyCoord(
            self.row["RAJ2000"],
            self.row["DEJ2000"],
            unit="deg",
            frame="icrs",
        )
        self.origin = self.row["z_class"]
        self.z_source = self.row["z_source"]
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
        spec_type = str(self.row["SpectrumType"])

        if "FHL" in str(self.row["Source_Name"]):
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
                f"Spectral model {spec_type!r} not implemented for FHL source"
            )

        if "FGL" in str(self.row["Source_Name"]):
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

        raise ValueError(
            f"Unknown Catalog for source name {self.name} for spectral model creation"
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


@dataclass(slots=True)
class OnOffSimulationInput:
    dataset_onoff: SpectrumDatasetOnOff
    npred_background: Any


class SourceAnalysis:
    """
    Concrete analysis of one source for one IRF node.

    This class combines one Source, one IRFNode, and observation-specific
    analysis settings. It loops over the redshift-specific spectral models owned
    by Source.
    """

    def __init__(
        self,
        source: Source,
        irf_node: IRFNode,
        output_table: QTable,
        *,
        offset: u.Quantity | None = None,  # Needs to match IRFs
        n_off_regions: int = 5,
        n_fake_datasets: int = 100,
    ):
        self.source = source
        self.irf_node = irf_node
        self.output_table = output_table
        self.offset = offset
        self.n_off_regions = n_off_regions
        self.n_fake_datasets = n_fake_datasets

        if offset is not None and offset.value > 0.0:
            self.pointing = self.source.position.directional_offset_by(
                position_angle=90 * u.deg,
                separation=self.offset,
            )
        else:
            self.pointing = self.source.position

        self.on_region = PointSkyRegion(self.source.position)

    def run(self) -> None:
        if self.source.spectral_models is None:
            raise ValueError("Source has no spectral model(s) to run analysis with")

        self.output_table["matched_node_zen"] = [self.irf_node.zen]
        self.output_table["matched_node_az"] = [self.irf_node.az]
        self.output_table["matched_node_delta_b"] = [self.irf_node.delta_b]
        self.output_table["matched_node_sin_delta"] = [self.irf_node.sin_delta]
        self.output_table["matched_node_cos_theta"] = [self.irf_node.cos_theta]

        observations = self._create_observations()
        for spectral_model in self.source.spectral_models:
            datasets_onoff = self.create_datasets_onoff(observations, spectral_model)
            sigma_results = self.run_fake_studies(datasets_onoff, spectral_model)
            obstime_results = self.estimate_obstime(sigma_results)
            self.append_results_to_output_table(obstime_results, spectral_model.name)

        self.output_table["status"] = [AnalysisStatus.SUCCESS]

    @staticmethod
    def _format_obstime(obstime: float) -> str:
        return f"{float(obstime):g}"

    def _create_observations(self) -> Observations:
        observations = {}
        for obs_id, obstime in enumerate(self.irf_node.obstimes, start=1):
            observations[obstime] = Observation.create(
                obs_id=obs_id,
                pointing=self.pointing,
                livetime=float(obstime) * u.h,
                irfs=self.irf_node.load_irfs(obstime),
                location=self.irf_node.location,
            )

        return observations

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

        rad_max = observation.rad_max.evaluate(
            energy=energy,
            offset=source_offset,
        )
        theta = rad_max.to(u.rad)
        solid_angle = 2.0 * np.pi * (1.0 - np.cos(theta.value)) * u.sr

        bkg_rate = observation.bkg.evaluate(
            energy=energy,
            offset=source_offset,
        )
        livetime = observation.observation_live_time_duration

        background_counts = (bkg_rate * energy_width * livetime * solid_angle).to_value(
            ""
        )
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
            observation,
            energy_axis_reco,
        )

        background_data = np.zeros(dataset.counts.data.shape, dtype=float)
        background_data[...] = background_counts.reshape(
            (energy_axis_reco.nbin,) + (1,) * (background_data.ndim - 1)
        )

        dataset.background = dataset.counts.copy(data=background_data)

    def create_spectrum_dataset_onoff(
        self,
        observation: Observation,
        spectral_model: SkyModel,
        obstime: float,
    ) -> OnOffSimulationInput:
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
            name=f"{self.source.name}_{spectral_model.name}_{obstime:g}h",
        )

        dataset_maker = SpectrumDatasetMaker(
            containment_correction=False, selection=["exposure", "edisp"]
        )
        safe_mask_maker = SafeMaskMaker(methods=["aeff-default"])

        dataset = dataset_maker.run(dataset_empty, observation)
        self._set_rad_max_background(dataset, observation, energy_axis_reco)
        dataset = safe_mask_maker.run(dataset, observation)
        dataset.models = spectral_model.copy()

        npred_background = dataset.npred_background()

        dataset_onoff = SpectrumDatasetOnOff.from_spectrum_dataset(
            dataset=dataset,
            acceptance=1,
            acceptance_off=self.n_off_regions,
        )

        return OnOffSimulationInput(
            dataset_onoff=dataset_onoff,
            npred_background=npred_background,
        )

    def create_datasets_onoff(
        self,
        observations: dict[float, Observations],
        spectral_model: SkyModel,
    ) -> dict[float, OnOffSimulationInput]:
        datasets_by_obstime: dict[float, OnOffSimulationInput] = {}
        for obstime, observation in observations.items():
            datasets_by_obstime[obstime] = self.create_spectrum_dataset_onoff(
                observation,
                spectral_model,
                obstime,
            )

        return datasets_by_obstime

    def run_fake_studies(
        self,
        datasets_onoff: dict[float, OnOffSimulationInput],
        spectral_model: SkyModel,
        *,
        n_fake_datasets: int | None = None,
    ) -> dict[float, dict[str, Any]]:
        results: dict[float, dict[str, Any]] = {}
        for obstime, dataset_on_off in datasets_onoff.items():
            fake_datasets, info_table, sigma = self.fake_data(
                dataset_on_off,
                spectral_model,
                n_fake_datasets=n_fake_datasets,
            )

            valid_sigma = sigma[np.isfinite(sigma)]
            if len(valid_sigma) == 0:
                sigma_mean = np.nan
                sigma_std = np.nan
            else:
                sigma_mean = float(np.mean(valid_sigma))
                sigma_std = float(np.std(valid_sigma))

            results[obstime] = {
                "dataset_onoff": dataset_on_off,
                "fake_datasets": fake_datasets,
                "info_table": info_table,
                "sigma": sigma,
                "sigma_mean": sigma_mean,
                "sigma_std": sigma_std,
            }

        return dict(sorted(results.items()))

    def fake_data(
        self,
        simulation_input: OnOffSimulationInput,
        spectral_model: SkyModel,
        *,
        n_fake_datasets: int | None = None,
    ) -> tuple[Datasets, QTable, np.ndarray]:
        n_fake = self.n_fake_datasets if n_fake_datasets is None else n_fake_datasets

        datasets = Datasets()
        for idx in range(n_fake):
            ds = simulation_input.dataset_onoff.copy(
                name=f"{simulation_input.dataset_onoff.name}_{idx}"
            )
            ds.models = spectral_model.copy()

            npred_background = simulation_input.npred_background.copy()
            data = npred_background.data

            invalid = ~np.isfinite(data) | (data < 0)
            if np.any(invalid):
                log.warning(
                    "%s: replacing %d invalid npred_background bins with 0",
                    ds.name,
                    int(np.count_nonzero(invalid)),
                )
                data[invalid] = 0.0

            ds.fake(random_state=idx, npred_background=npred_background)
            ds.meta_table["OBS_ID"] = [idx]
            datasets.append(ds)

        info_table = datasets.info_table()
        sigma = np.asarray(info_table["sqrt_ts"], dtype=float)

        return datasets, info_table, sigma

    def estimate_obstime(
        self,
        results: dict[float, dict[str, Any]],
    ) -> dict[str, Any]:
        sigma_by_obstime = {
            obstime: payload["sigma"] for obstime, payload in results.items()
        }

        obstime_5s, obstime_5s_std = self._fit_obstime_from_significance(
            sigma_by_obstime=sigma_by_obstime,
            sigma_target=self.output_table[0]["sigma_target"],
        )

        return {
            "obstime_5s": float(obstime_5s),
            "obstime_5s_std": float(obstime_5s_std),
            "results": results,
        }

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

        # S ~ sqrt(T)
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

    def append_results_to_output_table(
        self,
        obstime_results: dict[str, Any],
        model_label: str,
    ) -> None:
        label = ""
        if self.source.has_prior_redshift_scenarios:
            label = f"_{model_label}"

        self.output_table[f"obstime_5s{label}"] = obstime_results["obstime_5s"]
        self.output_table[f"obstime_5s_std{label}"] = obstime_results["obstime_5s_std"]
        for obstime, result in obstime_results["results"].items():
            obstime_label = self._format_obstime(obstime)
            self.output_table[f"sigma{label}_{obstime_label}h_mean"] = result[
                "sigma_mean"
            ]
            self.output_table[f"sigma{label}_{obstime_label}h_std"] = result[
                "sigma_std"
            ]

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
        """
        source_name = self.source.name

        if "FHL" in source_name:
            edges = [10, 20, 50, 150, 500, 2000] * u.GeV

        elif "FGL" in source_name:
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
            log.warning(
                "%s: cannot infer Fermi flux-point energy bins from source name",
                source_name,
            )
            return None

        if len(edges) != n_bins + 1:
            log.warning(
                "%s: expected %d flux-point edges for %d bins, got %d edges",
                source_name,
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

        catalog_label = "3FHL" if "FHL" in self.source.name else "4FGL"

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

        if "FGL" in self.source.name:
            e_lim = [5.0e-5, 1.0e3]
            if energy_bounds is None:
                energy_bounds = [5.0e-5, 100.0] * u.TeV
        elif "FHL" in self.source.name:
            e_lim = [5.0e-3, 1.0e3]
            if energy_bounds is None:
                energy_bounds = [5.0e-3, 100.0] * u.TeV
        else:
            e_lim = [5.0e-3, 1.0e3]
            if energy_bounds is None:
                energy_bounds = [5.0e-3, 100.0] * u.TeV

        fig, ax = plt.subplots()

        obstimes = self.irf_node.obstimes

        if len(obstimes) == 1:
            alphas = np.array([1.0])
        else:
            alphas = np.linspace(0.3, 1.0, len(obstimes))

        for obstime, alpha in zip(obstimes, alphas, strict=True):
            sens = self.irf_node.load_benchmark(obstime)
            energy_center = 0.5 * (sens["ENERG_LO"] + sens["ENERG_HI"])
            xerr = 0.5 * (sens["ENERG_HI"] - sens["ENERG_LO"])

            ax.errorbar(
                energy_center.flatten(),
                sens["ENERGY_FLUX_SENSITIVITY"].flatten(),
                xerr=xerr.flatten(),
                ls="",
                color=CTAO_COLORS["interstellar_indigo"],
                alpha=alpha,
                label=f"CTAO-N - {obstime:g}h",
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
    obstimes: list[float],
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
        output_table[f"obstime_5s_{label}"] = [np.nan]
        output_table[f"obstime_5s_std_{label}"] = [np.nan]

        for obstime in obstimes:
            obstime_label = f"{obstime:g}"
            output_table[f"sigma_{label}_{obstime_label}h_mean"] = [np.nan]
            output_table[f"sigma_{label}_{obstime_label}h_std"] = [np.nan]

    output_table["obstime_5s"] = [np.nan]
    output_table["obstime_5s_std"] = [np.nan]

    for obstime in obstimes:
        obstime_label = f"{obstime:g}"
        output_table[f"sigma_{obstime_label}h_mean"] = [np.nan]
        output_table[f"sigma_{obstime_label}h_std"] = [np.nan]

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
    output_table = create_output_table(
        source.table, irf_collection.obstimes, sigma_target=sigma_target
    )

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
