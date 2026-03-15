from astropy.coordinates import angular_separation
import astropy.units as u
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np

from ctapipe.reco.preprocessing import horizontal_to_telescope


def plot_theta_2(dl2_table):
    fig_theta_2, ax_theta_2 = plt.subplots(1, 1)

    alt = dl2_table["disp_alt"]
    az = dl2_table["disp_az"]
    lon, lat = horizontal_to_telescope(
        alt,
        az,
        dl2_table["subarray_pointing_lat"],
        dl2_table["subarray_pointing_lon"],
    )
    true_alt = dl2_table["true_alt"]
    true_az = dl2_table["true_az"]
    true_lon, true_lat = horizontal_to_telescope(
        true_alt,
        true_az,
        dl2_table["subarray_pointing_lat"],
        dl2_table["subarray_pointing_lon"],
    )
    off_angles = angular_separation(
        true_lon.to(u.rad),
        true_lat.to(u.rad),
        lon.to(u.rad),
        lat.to(u.rad),
    )
    theta_2 = (off_angles[np.isfinite(off_angles)]).to(u.deg) ** 2

    ax_theta_2.hist(theta_2, bins=50, histtype="step", color="blue")
    ax_theta_2.set_xlabel(r"$\theta^2$ / deg")
    ax_theta_2.set_ylabel("# of events")
    ax_theta_2.set_yscale("log")
    fig_theta_2.tight_layout()

    return fig_theta_2


def plot_reco_lon_lat(dl2_table):
    alt = dl2_table["disp_alt"]
    az = dl2_table["disp_az"]
    lon, lat = horizontal_to_telescope(
        alt,
        az,
        dl2_table["subarray_pointing_lat"],
        dl2_table["subarray_pointing_lon"],
    )

    lon_no_nans = lon[np.isfinite(lon)]
    lat_no_nans = lat[np.isfinite(lat)]

    fig_2d_hist, ax_2d_hist = plt.subplots(1, 1)
    max_abs_lon = np.max(np.abs(lon_no_nans.value))
    max_abs_lat = np.max(np.abs(lat_no_nans.value))
    max_abs = np.max([max_abs_lon, max_abs_lat])
    H, xedges, yedges = np.histogram2d(
        lon_no_nans.value,
        lat_no_nans.value,
        range=[
            [-max_abs, max_abs],
            [-max_abs, max_abs],
        ],
        bins=(50, 50),
    )

    im = ax_2d_hist.imshow(
        H.T,
        origin="lower",
        extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
        aspect="auto",
        cmap="viridis",
        norm=LogNorm(clip=True),
    )
    fig_2d_hist.colorbar(im, ax=ax_2d_hist, label="Counts")
    ax_2d_hist.set_xlabel("Longitude / deg")
    ax_2d_hist.set_ylabel("Latitude / deg")
    ax_2d_hist.set_title("Reconstructed Lon/Lat")
    ax_2d_hist.set_aspect("equal")
    fig_2d_hist.tight_layout()

    return fig_2d_hist


def stack_theta_hist(dl2_chunk):
    alt = dl2_chunk["disp_tel_alt"]
    az = dl2_chunk["disp_tel_az"]
    lon, lat = horizontal_to_telescope(
        alt,
        az,
        dl2_chunk["telescope_pointing_altitude"],
        dl2_chunk["telescope_pointing_azimuth"],
    )
    true_az = (dl2_chunk["true_az"],)
    true_alt = (dl2_chunk["true_alt"],)
    true_lon, true_lat = horizontal_to_telescope(
        true_alt,
        true_az,
        dl2_chunk["telescope_pointing_altitude"],
        dl2_chunk["telescope_pointing_azimuth"],
    )
    off_angles = angular_separation(
        true_lon.to(u.rad),
        true_lat.to(u.rad),
        lon.to(u.rad),
        lat.to(u.rad),
    )
    theta_2 = (off_angles[np.isfinite(off_angles)]).to(u.deg) ** 2
    hist, _ = np.histogram(theta_2, bins=100, range=u.Quantity([0, 30], u.deg**2))
    return hist


def stack_alt_az_hist(dl2_chunk):
    alt = dl2_chunk["disp_tel_alt"]
    az = dl2_chunk["disp_tel_az"]
    lon, lat = horizontal_to_telescope(
        alt,
        az,
        dl2_chunk["telescope_pointing_altitude"],
        dl2_chunk["telescope_pointing_azimuth"],
    )

    lon_no_nans = lon[np.isfinite(lon)]
    lat_no_nans = lat[np.isfinite(lat)]
    # az_mod = (az_no_nans + 180) % 360 - 180

    hist, _, _ = np.histogram2d(
        lon_no_nans.value, lat_no_nans.value, bins=(50, 50), range=[[-7, 7], [-7, 7]]
    )

    return hist


def plot_theta_2_tel(hist, tel_type):
    fig_theta_2, ax_theta_2 = plt.subplots(1, 1)

    bins = np.linspace(0, 30, 101)  # 50 bins → 51 edges

    ax_theta_2.step(bins[:-1], hist, where="mid", color="blue")
    ax_theta_2.set_title(rf"$\theta^2$  of {tel_type}")
    ax_theta_2.set_xlabel(r"$\theta^2$ / deg")
    ax_theta_2.set_ylabel("# of events")
    ax_theta_2.set_yscale("log")
    fig_theta_2.tight_layout()

    return fig_theta_2


def plot_reco_alt_az_tel(hist, tel_type):
    fig_2d_hist, ax_2d_hist = plt.subplots(1, 1)
    # az_edges = np.linspace(-15, 15, 51)  # 50 bins → 51 edges
    # alt_edges = np.linspace(55, 85, 51)
    h = ax_2d_hist.imshow(
        hist.T,
        origin="lower",
        extent=[-15, 15, 55, 85],
        cmap="viridis",
        norm=LogNorm(clip=True),
        aspect="equal",
    )

    ax_2d_hist.set_xlabel("Longitude / deg")
    ax_2d_hist.set_ylabel("Latitude / deg")
    ax_2d_hist.set_title(f"Reconstructed Lon/Lat for {tel_type}")
    ax_2d_hist.set_aspect("equal")
    fig_2d_hist.colorbar(h, ax=ax_2d_hist, label="Counts")
    fig_2d_hist.tight_layout()

    return fig_2d_hist
