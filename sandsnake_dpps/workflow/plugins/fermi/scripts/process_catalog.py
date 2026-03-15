import numpy as np
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_sun
from astropy.table import join, QTable
from astropy.time import Time
import astropy.units as u
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import argparse
import os
from tqdm import tqdm

astro_night = -18 * u.deg
# LaPalma geomag from MC's
la_palma = EarthLocation.of_site("Roque de los Muchachos")
prod_site_B_declination = -0.07534063 * u.rad
prod_site_B_inclination = 0.64818996 * u.rad
prod_site_B_toal = 38.77302551 * u.uT


def hist_alt(source_alt, observable_mask, alt_min, step_minutes, bin_width=5):
    """
    Plot a histogram of source altitudes during observable time (i.e., when the source is
    above the minimum altitude and the Sun is below -18 deg).

    The histogram shows observation time (in hours) spent per altitude bin.
    The function also computes and returns:
    - the median altitude,
    - the bin with the maximum accumulated observation time,
    - and the center of that bin.

    Parameters
    ----------
    source_alt : Quantity
        Altitude values of the source (in degrees) for each time step.
    observable_mask : ndarray of bool
        Boolean mask indicating when the source is observable (astronomical night and above threshold).
    alt_min : Quantity
        Minimum altitude threshold (e.g., 20 deg).
    step_minutes : int
        Time resolution of the time grid in minutes.
    bin_width : int, optional
        Width of the histogram bins in degrees (default is 5 deg).

    Returns
    -------
    fig : matplotlib.figure.Figure
        The generated matplotlib figure object.
    alt_median : float
        Median altitude (in degrees) during observable time, or NaN if no data.
    alt_max_hours : float
        Maximum total observation time (in hours) in any single altitude bin.
    alt_max_hours_bin : float
        Center of the altitude bin with the most observation time (in degrees).
    bin_width : float
        Width of the altitude histogram bins (in degrees).
    """
    observable_alt = source_alt[observable_mask]

    if len(observable_alt) != 0:
        alt_median = np.nanmedian(observable_alt)
    else:
        alt_median = np.nan * u.deg

    bin_edges = np.arange(alt_min.value, 90 + bin_width, bin_width)
    hist, _ = np.histogram(observable_alt.value, bins=bin_edges)

    # Convert to hours
    hist_hours = (hist * step_minutes).to(u.h)

    alt_max_hours = hist_hours.max()
    max_bin_index = np.argmax(hist_hours)
    alt_max_hours_bin = (bin_edges[max_bin_index] + bin_edges[max_bin_index + 1]) / 2

    fig, ax = plt.subplots()
    ax.bar(
        bin_edges[:-1],
        hist_hours.value,
        width=bin_width,
        align="edge",
        edgecolor="black",
    )
    ax.set_xlabel("Alt / deg")
    ax.set_ylabel("Obstime / h")
    ax.set_title("Altitude distribution during astronomical night")
    ax.grid(True)
    fig.tight_layout()

    if np.isfinite(alt_median):
        ax.axvline(
            alt_median.value,
            color="red",
            linestyle="--",
            label=f"Median Alt: {alt_median:.2f}",
        )
        ax.legend()

    return fig, alt_median, alt_max_hours, alt_max_hours_bin, bin_width


def hist_delta_B(delta_B, step_minutes, bin_width=10):
    """
    Plot a histogram of delta_B angles (angle between line-of-sight and geomagnetic field)
    weighted by observation time (in hours). Also computes and marks the median angle
    and identifies the bin with the highest accumulated observation time.

    Parameters
    ----------
    delta_B : Quantity
        Angles between line-of-sight (LOS) and geomagnetic field in degrees
        for all observable time steps.
    step_minutes : int
        Time resolution of the time grid in minutes.
    bin_width : int, optional
        Width of histogram bins in degrees (default is 10 deg).

    Returns
    -------
    fig : matplotlib.figure.Figure
        The generated matplotlib figure object.
    delta_B_median : float
        Median of delta_B angles (in degrees), or NaN if no data.
    alt_max_hours : float
        Maximum total observation time (in hours) in any single delta_B bin.
    alt_max_hours_bin : float
        Center of the bin with the most observation time (in degrees).
    bin_width : float
        Width of the delta_B histogram bins (in degrees).
    """
    if len(delta_B) != 0:
        delta_B_median = np.nanmedian(delta_B)
    else:
        delta_B_median = np.nan * u.deg

    bins = np.arange(0, 90 + bin_width, bin_width)
    hist_counts, bin_edges = np.histogram(delta_B.value, bins=bins)

    # Convert to hours
    hist_hours = (hist_counts * step_minutes).to(u.h)

    delta_B_max_hours = hist_hours.max()
    max_bin_index = np.argmax(hist_hours)
    delta_B_max_hours_bin = (
        bin_edges[max_bin_index] + bin_edges[max_bin_index + 1]
    ) / 2

    fig, ax = plt.subplots()
    ax.bar(
        bin_edges[:-1],
        hist_hours.value,
        width=bin_width,
        align="edge",
        edgecolor="black",
        alpha=0.8,
    )

    ax.set_xlabel("Angle between LOS and B field / deg")
    ax.set_ylabel("Observable time / h")
    ax.set_title(r"Distribution of $\Delta B$ angles")
    ax.grid(True)
    fig.tight_layout()
    if np.isfinite(delta_B_median):
        ax.axvline(
            delta_B_median.value,
            color="red",
            linestyle="--",
            label=f"median = {delta_B_median:.1f}",
        )

        ax.legend()

    return fig, delta_B_median, delta_B_max_hours, delta_B_max_hours_bin, bin_width


def calc_max_obstime(
    source_altaz,
    sun_altaz,
    alt_min,
    step_minutes,
):
    """
    Calculate the total observable time for a source above a given altitude
    and during astronomical night.

    Parameters
    ----------
    source_altaz : SkyCoord
        Source coordinates in AltAz frame over time.
    sun_altaz : SkyCoord
        Sun coordinates in AltAz frame over time.
    alt_min : Quantity
        Minimum required altitude (e.g. 20 deg).
    step_minutes : int
        Time resolution in minutes.

    Returns
    -------
    max_obs_time : Quantity
        Total observable time in hours.
    observable_mask : ndarray of bool
        Boolean mask indicating valid observation times.
    """
    is_night = sun_altaz.alt < astro_night
    is_visible = source_altaz.alt > alt_min
    observable_mask = is_night & is_visible

    max_obs_time = np.nan * u.h
    if observable_mask.any():
        max_obs_time = observable_mask.sum() * (step_minutes.to(u.h))
    return max_obs_time, observable_mask


def get_B_direction(dec, incl):
    """
    Compute the unit vector of the geomagnetic field in the local ENU system
    based on magnetic declination and inclination.

    Parameters
    ----------
    dec : Quantity
        Magnetic declination (angle east from geographic north) in radians.
    incl : Quantity
        Magnetic inclination (angle from horizontal, positive downward) in radians.

    Returns
    -------
    ndarray
        Normalized geomagnetic field direction vector [East, North, Up].
    """
    Bx = np.cos(incl) * np.sin(dec)
    By = np.cos(incl) * np.cos(dec)
    Bz = -np.sin(incl)

    return np.array([Bx, By, Bz]) / np.linalg.norm([Bx, By, Bz])


def calc_delta_B(source_altaz, observable_mask, B_decl, B_incl):
    """
    Compute the angle (delta_B) between line-of-sight to the source and
    the local geomagnetic field vector during observable times.

    Parameters
    ----------
    source_altaz : SkyCoord
        AltAz coordinates of the source over time.
    observable_mask : ndarray of bool
        Boolean mask indicating observable time steps.
    B_decl : Quantity
        Magnetic declination (radians).
    B_incl : Quantity
        Magnetic inclination (radians).

    Returns
    -------
    delta_B : Quantity
        Array of angles in degrees between LOS and geomagnetic field,
        centered around 90° (i.e., angle relative to perpendicularity).
    """
    observable_altaz = source_altaz[observable_mask]
    B_vec = get_B_direction(B_decl, B_incl)

    # LOS vectors, shape (3, N)
    los_xyz = observable_altaz.cartesian.xyz.value
    cos_theta = np.dot(los_xyz.T, B_vec)

    delta_B = np.arccos(cos_theta) * u.rad
    delta_B = np.abs(90 * u.deg - delta_B.to(u.deg))
    return delta_B


def trim_source_names(table, columns):
    for col in columns:
        if col == "ASSOC_FHL":
            table[col] = [
                str(sn.strip()).replace(" ", "_") if sn is not np.ma.masked else ""
                for sn in table[col]
            ]
        else:
            table[col] = [str(sn.strip()).replace(" ", "_") for sn in table[col]]
    return table, table[columns]


def main(catalog, fgl_path, lac_path, fhl_path, outdir):
    """
    Process a Fermi-LAT catalog and compute visibility-related parameters for each source
    assuming observation from CTAO-North (La Palma) for a
    specific time interval (~ Obscycle 26/27).

    For each source:
    - Compute total observable time (above 20° altitude during astronomical night)
    - Compute median altitude during observable time
    - Compute angle between line-of-sight and local geomagnetic field (delta_B)
    - Compute median delta_B angle
    - Generate altitude and delta_B histograms as PDF
    - Save per-source results in individual ECSV files
    """
    outpath = Path(outdir)
    outpath.mkdir(parents=True, exist_ok=True)

    for cat_path in [fgl_path, lac_path, fhl_path]:
        if not os.path.exists(cat_path):
            raise FileNotFoundError(f"Catalog file not found: {cat_path}")

    catalog_dict = {"FGL": "4FGL_DR4", "LAC": "4LAC_DR3", "FHL": "3FHL_DR3"}
    if catalog not in list(catalog_dict.values()):
        raise ValueError(
            f"Not a valid catalog: {catalog}. "
            f"Valid ones are {list(catalog_dict.values())}."
        )

    columns_to_append = [
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
        "is_LAC_Source",  # Added column
    ]
    lac_table = QTable.read(lac_path, hdu=1)
    fgl_table = QTable.read(fgl_path, hdu=1)
    lac_table, _ = trim_source_names(lac_table, ["Source_Name"])
    fgl_table, src_names_fgl_assoc_fhl = trim_source_names(
        fgl_table, ["Source_Name", "ASSOC_FHL"]
    )

    if "FGL" in catalog:
        # Merge some LAC parameters to FGL catalog
        catalog_table = QTable.read(fgl_path, hdu=1)
        lac_table["is_LAC_Source"] = True
        lac_subtable = lac_table[columns_to_append]
        catalog_table, _ = trim_source_names(
            catalog_table, ["Source_Name", "ASSOC_FHL"]
        )
        catalog_table = join(
            catalog_table,
            lac_subtable,
            keys="Source_Name",
            join_type="left",
            table_names=["FGL", "LAC"],
        )
        for col in columns_to_append:
            if col not in ["Source_Name", "is_LAC_Source"]:
                catalog_table[col] = catalog_table[col].filled(np.nan)
            elif col == "is_LAC_Source":
                catalog_table[col] = catalog_table[col].filled(False)

    elif "FHL" in catalog:
        catalog_table = QTable.read(fhl_path, hdu=1)
        catalog_table, _ = trim_source_names(catalog_table, ["Source_Name"])

    step_minutes = 20 * u.min
    start = "2024-04-01"
    end = "2025-04-01"
    alt_min = 20 * u.deg
    time_grid = Time(
        np.arange(Time(start).jd, Time(end).jd, step_minutes.to(u.d).value), format="jd"
    )
    altaz_frame = AltAz(obstime=time_grid, location=la_palma)
    sun_altaz = get_sun(time_grid).transform_to(altaz_frame)

    for source in tqdm(catalog_table[:10], desc="Processing sources", unit="src"):
        source_name = source["Source_Name"]
        source_table = QTable(rows=[source], names=catalog_table.colnames)
        ra = source["RAJ2000"]
        dec = source["DEJ2000"]

        source_altaz = SkyCoord(ra=ra, dec=dec).transform_to(altaz_frame)

        max_obs_time, observable_mask = calc_max_obstime(
            source_altaz, sun_altaz, alt_min, step_minutes
        )

        if max_obs_time is not np.nan:
            delta_B = calc_delta_B(
                source_altaz,
                observable_mask,
                prod_site_B_declination,
                prod_site_B_inclination,
            )
            (
                fig_hist_alt,
                alt_median,
                alt_max_hours,
                alt_max_hours_bin,
                alt_bin_width,
            ) = hist_alt(
                source_altaz.alt,
                observable_mask,
                alt_min,
                step_minutes,
            )
            (
                fig_hist_delta_B,
                delta_B_median,
                delta_B_max_hours,
                delta_B_max_hours_bin,
                delta_B_bin_width,
            ) = hist_delta_B(delta_B, step_minutes)

            # check dir
            src_dir = outpath / source_name
            src_dir.mkdir(parents=True, exist_ok=True)

            with PdfPages(
                outpath / f"{source_name}/alt_delta_B_hist_{source_name}.pdf"
            ) as pdf:
                pdf.savefig(fig_hist_alt)
                pdf.savefig(fig_hist_delta_B)
        else:
            alt_median = delta_B_median = np.nan * u.deg
            alt_max_hours = delta_B_max_hours = np.nan * u.h
            alt_max_hours_bin = alt_bin_width = delta_B_max_hours_bin = (
                delta_B_bin_width
            ) = np.nan

        source_table["max_obstime"] = max_obs_time
        source_table["alt_median"] = alt_median
        source_table["alt_max_hours"] = alt_max_hours
        source_table["alt_max_hours_bin"] = alt_max_hours_bin
        source_table["alt_bin_width"] = alt_bin_width
        source_table["delta_B_median"] = delta_B_median
        source_table["delta_B_max_hours"] = delta_B_max_hours
        source_table["delta_B_max_hours_bin"] = delta_B_max_hours_bin
        source_table["delta_B_bin_width"] = delta_B_bin_width

        filename = outpath / f"{source_name}/{source_name}.ecsv"
        os.makedirs(filename.parent, exist_ok=True)
        source_table.write(filename, format="ascii.ecsv", overwrite=True)

        # create symlinks
        if "FGL" in catalog:
            if bool(source["ASSOC_FHL"]):
                dest_dir = outpath.parent / catalog_dict["FHL"] / source["ASSOC_FHL"]
                link_dir = outpath / source_name / source["ASSOC_FHL"]
                os.symlink(dest_dir, link_dir)

        if "FHL" in catalog:
            src_names_fhl = src_names_fgl_assoc_fhl["ASSOC_FHL"].tolist()
            src_names_fgl = src_names_fgl_assoc_fhl["Source_Name"].tolist()
            if source_name in src_names_fhl:
                src_name_fgl = src_names_fgl[src_names_fhl.index(source_name)]
                dest_dir = outpath.parent / catalog_dict["FGL"] / src_name_fgl
                link_dir = outpath / source_name / src_name_fgl
                os.symlink(dest_dir, link_dir)


if __name__ == "__main__":
    # TODO: 4FHL recently announced but not yet published.
    parser = argparse.ArgumentParser()
    parser.add_argument("--catalog", required=True)
    parser.add_argument("--fgl", required=True)
    parser.add_argument("--lac", required=True)
    parser.add_argument("--fhl", required=True)
    parser.add_argument("-o", "--outdir", required=True)
    args = parser.parse_args()

    main(args.catalog, args.fgl, args.lac, args.fhl, args.outdir)
