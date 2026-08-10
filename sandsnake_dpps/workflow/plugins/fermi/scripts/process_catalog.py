import argparse
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Self

import astropy.units as u
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
from astropy.coordinates import AltAz, EarthLocation, SkyCoord, get_sun
from astropy.table import QTable, join
from astropy.time import Time
from tqdm import tqdm

try:
    from plugins.fermi.scripts.catalog_priors import (
        RedshiftPriorConfig,
        append_redshift_prior_columns,
        build_redshift_prior_table,
    )
except ImportError:  # allows direct execution from plugins/fermi/scripts
    from catalog_priors import (  # type: ignore[no-redef]
        RedshiftPriorConfig,
        append_redshift_prior_columns,
        build_redshift_prior_table,
    )


CATALOG_NAMES = {
    "FGL": "4FGL_DR4",
    "LAC": "4LAC_DR3",
    "FHL": "3FHL_DR3",
}

LAC_COLUMNS_TO_APPEND = [
    "Source_Name",  # join key
    "VLBI_Counterpart",
    "Redshift",
    "SED_class",
    "HE_EPeak",
    "Unc_HE_EPeak",
    "HE_nuFnuPeak",
    "Unc_HE_nuFnuPeak",
    "nu_syn",
    "nuFnu_syn",
    "Highest_energy",
    "LAC_sample",
    "is_LAC_Source",
]


# TODO: Make prod site configurable
@dataclass(slots=True)
class VisibilityConfig:
    start: str = "2027-04-01"
    end: str = "2028-04-01"
    step_minutes: u.Quantity = 20 * u.min
    alt_min: u.Quantity = 60 * u.deg
    site_name: str = "Roque de los Muchachos"
    astro_night: u.Quantity = -18 * u.deg
    # La Palma geomagnetic field values from the MC production site configuration.
    prod_site_B_declination: u.Quantity = -0.07534063 * u.rad
    prod_site_B_inclination: u.Quantity = 0.64818996 * u.rad
    prod_site_B_total: u.Quantity = 38.77302551 * u.uT


@dataclass(slots=True)
class VisibilityGrid:
    config: VisibilityConfig
    location: EarthLocation
    time_grid: Time
    altaz_frame: AltAz
    sun_altaz: SkyCoord
    is_night: np.ndarray
    b_vec: np.ndarray

    @classmethod
    def create(cls, config: VisibilityConfig) -> Self:
        location = EarthLocation.of_site(config.site_name)
        time_grid = Time(
            np.arange(
                Time(config.start).jd,
                Time(config.end).jd,
                config.step_minutes.to_value(u.d),
            ),
            format="jd",
        )
        altaz_frame = AltAz(obstime=time_grid, location=location)
        sun_altaz = get_sun(time_grid).transform_to(altaz_frame)
        is_night = np.asarray(sun_altaz.alt < config.astro_night, dtype=bool)
        b_vec = get_B_direction(
            config.prod_site_B_declination, config.prod_site_B_inclination
        )

        return cls(
            config=config,
            location=location,
            time_grid=time_grid,
            altaz_frame=altaz_frame,
            sun_altaz=sun_altaz,
            is_night=is_night,
            b_vec=b_vec,
        )


@dataclass(slots=True)
class VisibilityResult:
    max_obstime: u.Quantity
    cos_theta_mean: float
    alt_max_hours: u.Quantity
    alt_max_hours_bin: float
    alt_bin_width: float
    sin_delta_mean: float
    delta_max_hours: u.Quantity
    delta_max_hours_bin: float
    delta_bin_width: float
    fig_hist_alt: Any | None = None
    fig_hist_delta_B: Any | None = None


def empty_visibility_result() -> VisibilityResult:
    return VisibilityResult(
        max_obstime=np.nan * u.h,
        cos_theta_mean=np.nan,
        alt_max_hours=np.nan * u.h,
        alt_max_hours_bin=np.nan,
        alt_bin_width=np.nan,
        sin_delta_mean=np.nan,
        delta_max_hours=np.nan * u.h,
        delta_max_hours_bin=np.nan,
        delta_bin_width=np.nan,
        fig_hist_alt=None,
        fig_hist_delta_B=None,
    )


class SourceVisibilityCalculator:
    def __init__(self, grid: VisibilityGrid, *, write_plots: bool = False):
        self.grid = grid
        self.write_plots = write_plots

    def calculate(self, row) -> VisibilityResult:
        source_position = SkyCoord(
            ra=row["RAJ2000"],
            dec=row["DEJ2000"],
            unit="deg",
            frame="icrs",
        )
        source_altaz = source_position.transform_to(self.grid.altaz_frame)

        max_obstime, observable_mask = calc_max_obstime(
            source_altaz=source_altaz,
            is_night=self.grid.is_night,
            alt_min=self.grid.config.alt_min,
            step_minutes=self.grid.config.step_minutes,
        )

        if not observable_mask.any():
            return empty_visibility_result()

        delta_B = calc_delta_B(
            source_altaz=source_altaz,
            observable_mask=observable_mask,
            B_decl=self.grid.config.prod_site_B_declination,
            B_incl=self.grid.config.prod_site_B_inclination,
            B_vec=self.grid.b_vec,
        )

        alt_summary = hist_alt(
            source_alt=source_altaz.alt,
            observable_mask=observable_mask,
            alt_min=self.grid.config.alt_min,
            step_minutes=self.grid.config.step_minutes,
            make_plot=self.write_plots,
        )
        delta_summary = hist_delta_B(
            delta_B=delta_B,
            step_minutes=self.grid.config.step_minutes,
            make_plot=self.write_plots,
        )

        return VisibilityResult(
            max_obstime=max_obstime,
            cos_theta_mean=alt_summary.mean_value,
            alt_max_hours=alt_summary.max_hours,
            alt_max_hours_bin=alt_summary.max_hours_bin,
            alt_bin_width=alt_summary.bin_width,
            sin_delta_mean=delta_summary.mean_value,
            delta_max_hours=delta_summary.max_hours,
            delta_max_hours_bin=delta_summary.max_hours_bin,
            delta_bin_width=delta_summary.bin_width,
            fig_hist_alt=alt_summary.fig,
            fig_hist_delta_B=delta_summary.fig,
        )


@dataclass(slots=True)
class HistogramSummary:
    fig: Any | None
    mean_value: float
    max_hours: u.Quantity
    max_hours_bin: float
    bin_width: float


def hist_alt(
    source_alt: u.Quantity,
    observable_mask: np.ndarray,
    alt_min: u.Quantity,
    step_minutes: u.Quantity,
    bin_width: float = 5,
    *,
    make_plot: bool = True,
) -> HistogramSummary:
    """
    Summarize source altitude during observable time.

    The mean value is mean(cos(theta)), where theta = 90 deg - altitude.
    """
    observable_alt = source_alt[observable_mask]

    if len(observable_alt) != 0:
        zen = (90.0 * u.deg - observable_alt).to_value(u.rad)
        cos_theta_mean = float(np.mean(np.cos(zen)))
    else:
        cos_theta_mean = np.nan

    bin_edges = np.arange(alt_min.to_value(u.deg), 90 + bin_width, bin_width)
    hist, _ = np.histogram(observable_alt.to_value(u.deg), bins=bin_edges)
    hist_hours = (hist * step_minutes).to(u.h)

    if len(hist_hours) == 0 or hist_hours.sum() == 0 * u.h:
        alt_max_hours = np.nan * u.h
        alt_max_hours_bin = np.nan
    else:
        max_bin_index = int(np.argmax(hist_hours))
        alt_max_hours = hist_hours[max_bin_index]
        alt_max_hours_bin = float(
            0.5 * (bin_edges[max_bin_index] + bin_edges[max_bin_index + 1])
        )

    fig = None
    if make_plot:
        fig, ax = plt.subplots()
        ax.bar(
            bin_edges[:-1],
            hist_hours.to_value(u.h),
            width=bin_width,
            align="edge",
            edgecolor="black",
        )
        ax.set_xlabel("Alt / deg")
        ax.set_ylabel("Obstime / h")
        ax.set_title("Altitude distribution during astronomical night")
        ax.grid(True)

    return HistogramSummary(
        fig=fig,
        mean_value=cos_theta_mean,
        max_hours=alt_max_hours,
        max_hours_bin=alt_max_hours_bin,
        bin_width=float(bin_width),
    )


def hist_delta_B(
    delta_B: u.Quantity,
    step_minutes: u.Quantity,
    bin_width: float = 10,
    *,
    make_plot: bool = True,
) -> HistogramSummary:
    """
    Summarize delta_B angles during observable time.

    The mean value is mean(sin(delta_B)).
    """
    if len(delta_B) != 0:
        sin_delta_mean = float(np.mean(np.sin(delta_B.to_value(u.rad))))
    else:
        sin_delta_mean = np.nan

    bin_edges = np.arange(0, 90 + bin_width, bin_width)
    hist_counts, _ = np.histogram(delta_B.to_value(u.deg), bins=bin_edges)
    hist_hours = (hist_counts * step_minutes).to(u.h)

    if len(hist_hours) == 0 or hist_hours.sum() == 0 * u.h:
        delta_B_max_hours = np.nan * u.h
        delta_B_max_hours_bin = np.nan
    else:
        max_bin_index = int(np.argmax(hist_hours))
        delta_B_max_hours = hist_hours[max_bin_index]
        delta_B_max_hours_bin = float(
            0.5 * (bin_edges[max_bin_index] + bin_edges[max_bin_index + 1])
        )

    fig = None
    if make_plot:
        fig, ax = plt.subplots()
        ax.bar(
            bin_edges[:-1],
            hist_hours.to_value(u.h),
            width=bin_width,
            align="edge",
            edgecolor="black",
            alpha=0.8,
        )
        ax.set_xlabel("Angle between LOS and B field / deg")
        ax.set_ylabel("Observable time / h")
        ax.set_title(r"Distribution of $\Delta B$ angles")
        ax.grid(True)

    return HistogramSummary(
        fig=fig,
        mean_value=sin_delta_mean,
        max_hours=delta_B_max_hours,
        max_hours_bin=delta_B_max_hours_bin,
        bin_width=float(bin_width),
    )


def normalize_missing_string(value) -> str:
    if value is np.ma.masked:
        return ""
    text = str(value).strip()
    if text.lower() in {"", "--", "nan", "none", "masked"}:
        return ""
    return text


def clean_source_name(value) -> str:
    return normalize_missing_string(value).replace(" ", "_")


def trim_source_names(table: QTable, columns: list[str]) -> tuple[QTable, QTable]:
    for col in columns:
        if col not in table.colnames:
            continue
        table[col] = [clean_source_name(value) for value in table[col]]
    return table, table[[col for col in columns if col in table.colnames]]


def ensure_column(table: QTable, column: str, default_value) -> None:
    if column not in table.colnames:
        table[column] = [default_value] * len(table)


def ensure_redshift_prior_input_columns(table: QTable) -> None:
    """Ensure catalog_priors.py can run on both 4FGL-like and 3FHL-like tables."""
    if "CLASS1" not in table.colnames:
        if "CLASS" in table.colnames:
            table["CLASS1"] = table["CLASS"]
        elif "Source_Class" in table.colnames:
            table["CLASS1"] = table["Source_Class"]
        else:
            table["CLASS1"] = [""] * len(table)

    ensure_column(table, "SED_class", "")
    ensure_column(table, "Redshift", np.nan)


def calc_max_obstime(
    source_altaz: SkyCoord,
    is_night: np.ndarray,
    alt_min: u.Quantity,
    step_minutes: u.Quantity,
) -> tuple[u.Quantity, np.ndarray]:
    """Calculate total observable time above alt_min during astronomical night."""
    is_visible = np.asarray(source_altaz.alt > alt_min, dtype=bool)
    observable_mask = is_night & is_visible

    if observable_mask.any():
        return observable_mask.sum() * step_minutes.to(u.h), observable_mask

    return np.nan * u.h, observable_mask


def get_B_direction(dec: u.Quantity, incl: u.Quantity) -> np.ndarray:
    """Compute normalized geomagnetic field direction in local ENU coordinates."""
    Bx = np.cos(incl) * np.sin(dec)
    By = np.cos(incl) * np.cos(dec)
    Bz = -np.sin(incl)
    B = np.array([Bx, By, Bz], dtype=float)
    return B / np.linalg.norm(B)


def calc_delta_B(
    source_altaz: SkyCoord,
    observable_mask: np.ndarray,
    B_decl: u.Quantity,
    B_incl: u.Quantity,
    B_vec: np.ndarray | None = None,
) -> u.Quantity:
    """
    Compute delta_B for observable time steps.

    delta_B is abs(90 deg - angle(LOS, B)), matching the existing convention in
    the Fermi workflow.
    """
    observable_altaz = source_altaz[observable_mask]
    if len(observable_altaz) == 0:
        return np.array([], dtype=float) * u.deg

    if B_vec is None:
        B_vec = get_B_direction(B_decl, B_incl)

    los_xyz = np.asarray(observable_altaz.cartesian.xyz.value, dtype=float)
    cos_angle = np.clip(np.dot(los_xyz.T, B_vec), -1.0, 1.0)
    angle = np.arccos(cos_angle) * u.rad
    return np.abs(90.0 * u.deg - angle.to(u.deg))


def read_and_prepare_catalogs(
    catalog: str,
    fgl_path: str | Path,
    lac_path: str | Path,
    fhl_path: str | Path,
) -> tuple[QTable, QTable]:
    for cat_path in [fgl_path, lac_path, fhl_path]:
        if not Path(cat_path).exists():
            raise FileNotFoundError(f"Catalog file not found: {cat_path}")

    if catalog not in {CATALOG_NAMES["FGL"], CATALOG_NAMES["FHL"]}:
        raise ValueError(
            f"Not a valid catalog: {catalog}. "
            f"Valid ones are {[CATALOG_NAMES['FGL'], CATALOG_NAMES['FHL']]}"
        )

    lac_table = QTable.read(lac_path, hdu=1)
    fgl_table = QTable.read(fgl_path, hdu=1)

    lac_table, _ = trim_source_names(lac_table, ["Source_Name"])
    fgl_table, src_names_fgl_assoc_fhl = trim_source_names(
        fgl_table,
        ["Source_Name", "ASSOC_FHL"],
    )

    if catalog == CATALOG_NAMES["FGL"]:
        catalog_table = fgl_table.copy()

        lac_table["is_LAC_Source"] = True
        lac_subtable = lac_table[
            [col for col in LAC_COLUMNS_TO_APPEND if col in lac_table.colnames]
        ]

        catalog_table = join(
            catalog_table,
            lac_subtable,
            keys="Source_Name",
            join_type="left",
            table_names=["FGL", "LAC"],
        )

        for col in LAC_COLUMNS_TO_APPEND:
            if col not in catalog_table.colnames:
                continue
            if col == "Source_Name":
                continue
            if col == "is_LAC_Source":
                catalog_table[col] = catalog_table[col].filled(False)
            else:
                catalog_table[col] = catalog_table[col].filled(np.nan)

    else:
        catalog_table = QTable.read(fhl_path, hdu=1)
        catalog_table, _ = trim_source_names(catalog_table, ["Source_Name"])

    ensure_redshift_prior_input_columns(catalog_table)
    return catalog_table[:], src_names_fgl_assoc_fhl[:]


def add_redshift_priors(
    catalog_table: QTable,
    outpath: Path,
    *,
    lower_quantile: float = 0.16,
    upper_quantile: float = 0.84,
    min_sources_per_group: int = 50,
) -> QTable:
    ensure_redshift_prior_input_columns(catalog_table)

    prior_config = RedshiftPriorConfig(
        class_column="CLASS1",
        sed_column="SED_class",
        redshift_column="Redshift",
        lower_quantile=lower_quantile,
        upper_quantile=upper_quantile,
        min_sources_per_group=min_sources_per_group,
    )

    redshift_prior_table = build_redshift_prior_table(
        catalog_table,
        prior_config,
    )

    redshift_prior_table.write(
        outpath / "redshift_priors.ecsv",
        format="ascii.ecsv",
        overwrite=True,
    )

    return append_redshift_prior_columns(
        catalog_table,
        redshift_prior_table,
        prior_config,
    )


def append_visibility_columns(source_table: QTable, result: VisibilityResult) -> QTable:
    source_table["max_obstime"] = [result.max_obstime]
    source_table["cos_theta_mean"] = [result.cos_theta_mean]
    source_table["alt_max_hours"] = [result.alt_max_hours]
    source_table["alt_max_hours_bin"] = [result.alt_max_hours_bin]
    source_table["alt_bin_width"] = [result.alt_bin_width]
    source_table["sin_delta_mean"] = [result.sin_delta_mean]
    source_table["delta_max_hours"] = [result.delta_max_hours]
    source_table["delta_max_hours_bin"] = [result.delta_max_hours_bin]
    source_table["delta_bin_width"] = [result.delta_bin_width]
    return source_table


def write_histograms(
    source_dir: Path, source_name: str, result: VisibilityResult
) -> None:
    if result.fig_hist_alt is None or result.fig_hist_delta_B is None:
        return

    with PdfPages(source_dir / f"alt_delta_hist_{source_name}.pdf") as pdf:
        pdf.savefig(result.fig_hist_alt)
        pdf.savefig(result.fig_hist_delta_B)

    plt.close(result.fig_hist_alt)
    plt.close(result.fig_hist_delta_B)


def symlink_force(target: Path, link: Path) -> None:
    link.parent.mkdir(parents=True, exist_ok=True)
    if os.path.lexists(link):
        return
    link.symlink_to(target, target_is_directory=True)


def create_catalog_symlinks(
    *,
    catalog: str,
    source_row,
    source_name: str,
    outpath: Path,
    src_names_fgl_assoc_fhl: QTable,
) -> None:
    if catalog == CATALOG_NAMES["FGL"]:
        if "ASSOC_FHL" not in source_row.colnames:
            return
        assoc_fhl = clean_source_name(source_row["ASSOC_FHL"])
        if assoc_fhl:
            dest_dir = outpath.parent / CATALOG_NAMES["FHL"] / assoc_fhl
            link_dir = outpath / source_name / assoc_fhl
            symlink_force(dest_dir, link_dir)

    if catalog == CATALOG_NAMES["FHL"]:
        if "ASSOC_FHL" not in src_names_fgl_assoc_fhl.colnames:
            return
        src_names_fhl = [
            clean_source_name(value) for value in src_names_fgl_assoc_fhl["ASSOC_FHL"]
        ]
        src_names_fgl = [
            clean_source_name(value) for value in src_names_fgl_assoc_fhl["Source_Name"]
        ]
        if source_name in src_names_fhl:
            idx = src_names_fhl.index(source_name)
            src_name_fgl = src_names_fgl[idx]
            if src_name_fgl:
                dest_dir = outpath.parent / CATALOG_NAMES["FGL"] / src_name_fgl
                link_dir = outpath / source_name / src_name_fgl
                symlink_force(dest_dir, link_dir)


def process_sources(
    *,
    catalog: str,
    catalog_table: QTable,
    src_names_fgl_assoc_fhl: QTable,
    outpath: Path,
    calculator: SourceVisibilityCalculator,
    write_plots: bool,
) -> None:
    for source in tqdm(
        catalog_table,
        desc=f"Processing {catalog} sources",
        unit="src",
    ):
        source_name = clean_source_name(source["Source_Name"])
        if not source_name:
            raise ValueError("Encountered source without Source_Name:")

        source_table = QTable(rows=[source], names=catalog_table.colnames)
        result = calculator.calculate(source)
        source_table = append_visibility_columns(source_table, result)

        source_dir = outpath / source_name
        source_dir.mkdir(parents=True, exist_ok=True)

        if write_plots:
            write_histograms(source_dir, source_name, result)

        source_table.write(
            source_dir / f"{source_name}.ecsv",
            format="ascii.ecsv",
            overwrite=True,
        )

        create_catalog_symlinks(
            catalog=catalog,
            source_row=source,
            source_name=source_name,
            outpath=outpath,
            src_names_fgl_assoc_fhl=src_names_fgl_assoc_fhl,
        )
    return


def main(
    catalog: str,
    fgl_path: str | Path,
    lac_path: str | Path,
    fhl_path: str | Path,
    outdir: str | Path,
    *,
    start: str = "2027-04-01",
    end: str = "2028-04-01",
    step_minutes: float = 20.0,
    alt_min: float = 60.0,
    write_plots: bool = False,
    lower_quantile: float = 0.16,
    upper_quantile: float = 0.84,
    min_sources_per_group: int = 50,
) -> None:
    """
    Process a Fermi-LAT catalog and compute per-source catalog products.

    This step is intentionally catalog-side only:
    - read and merge catalog metadata,
    - derive redshift priors from sources with measured redshift,
    - precompute visibility summaries for CTAO-North,
    - write one enriched ECSV file per source.

    The later significance script should only consume these enriched source files.
    """
    outpath = Path(outdir)
    outpath.mkdir(parents=True, exist_ok=True)

    catalog_table, src_names_fgl_assoc_fhl = read_and_prepare_catalogs(
        catalog=catalog,
        fgl_path=fgl_path,
        lac_path=lac_path,
        fhl_path=fhl_path,
    )

    catalog_table = add_redshift_priors(
        catalog_table,
        outpath,
        lower_quantile=lower_quantile,
        upper_quantile=upper_quantile,
        min_sources_per_group=min_sources_per_group,
    )

    config = VisibilityConfig(
        start=start,
        end=end,
        step_minutes=step_minutes * u.min,
        alt_min=alt_min * u.deg,
    )
    grid = VisibilityGrid.create(config)
    calculator = SourceVisibilityCalculator(grid, write_plots=write_plots)

    process_sources(
        catalog=catalog,
        catalog_table=catalog_table,
        src_names_fgl_assoc_fhl=src_names_fgl_assoc_fhl,
        outpath=outpath,
        calculator=calculator,
        write_plots=write_plots,
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Prepare per-source Fermi catalog files for the CTAO significance workflow."
    )
    parser.add_argument("--catalog", required=True, choices=["4FGL_DR4", "3FHL_DR3"])
    parser.add_argument("--fgl", required=True)
    parser.add_argument("--lac", required=True)
    parser.add_argument("--fhl", required=True)
    parser.add_argument("-o", "--outdir", required=True)

    parser.add_argument("--start", default="2027-04-01")
    parser.add_argument("--end", default="2028-04-01")
    parser.add_argument("--step-minutes", type=float, default=20.0)
    parser.add_argument("--alt-min", type=float, default=60.0)
    parser.add_argument("--write-plots", action="store_true")

    parser.add_argument("--redshift-lower-quantile", type=float, default=0.16)
    parser.add_argument("--redshift-upper-quantile", type=float, default=0.84)
    parser.add_argument("--redshift-min-sources-per-group", type=int, default=50)

    return parser.parse_args()


def main_from_snakemake(snakemake) -> None:
    main(
        catalog=snakemake.wildcards.catalog,
        fgl_path=snakemake.input.fgl,
        lac_path=snakemake.input.lac,
        fhl_path=snakemake.input.fhl,
        outdir=snakemake.params.outdir,
        start=snakemake.params.start,
        end=snakemake.params.end,
        step_minutes=snakemake.params.step_minutes,
        alt_min=snakemake.params.alt_min,
        write_plots=snakemake.params.write_plots,
        lower_quantile=snakemake.params.redshift_lower_quantile,
        upper_quantile=snakemake.params.redshift_upper_quantile,
        min_sources_per_group=snakemake.params.redshift_min_sources_per_group,
    )


def main_from_args(args: argparse.Namespace) -> None:
    main(
        catalog=args.catalog,
        fgl_path=args.fgl,
        lac_path=args.lac,
        fhl_path=args.fhl,
        outdir=args.outdir,
        start=args.start,
        end=args.end,
        step_minutes=args.step_minutes,
        alt_min=args.alt_min,
        write_plots=args.write_plots,
        lower_quantile=args.redshift_lower_quantile,
        upper_quantile=args.redshift_upper_quantile,
        min_sources_per_group=args.redshift_min_sources_per_group,
    )


if "snakemake" in globals():
    main_from_snakemake(snakemake)  # noqa: F821
elif __name__ == "__main__":
    args = parse_args()
    main_from_args(args)
