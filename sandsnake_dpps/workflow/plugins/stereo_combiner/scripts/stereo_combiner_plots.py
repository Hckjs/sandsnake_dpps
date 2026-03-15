from astropy.coordinates import angular_separation
from astropy.table import QTable
import astropy.units as u
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
import math

from gammapy.irf import EffectiveAreaTable2D

from ctapipe.reco.preprocessing import horizontal_to_telescope
from ..mc.irf_plots import add_sensitivity_comparisons

e_lim = [5.0e-3, 5.0e2]


def stack_theta_hist(dl2_chunk):
    alt = dl2_chunk["disp_alt"]
    az = dl2_chunk["disp_az"]
    lon, lat = horizontal_to_telescope(
        alt,
        az,
        dl2_chunk["subarray_pointing_lat"],
        dl2_chunk["subarray_pointing_lon"],
    )
    true_az = dl2_chunk["true_az"]
    true_alt = dl2_chunk["true_alt"]
    true_lon, true_lat = horizontal_to_telescope(
        true_alt,
        true_az,
        dl2_chunk["subarray_pointing_lat"],
        dl2_chunk["subarray_pointing_lon"],
    )
    off_angles = angular_separation(
        true_lon.to(u.rad),
        true_lat.to(u.rad),
        lon.to(u.rad),
        lat.to(u.rad),
    )
    theta_2 = (off_angles[np.isfinite(off_angles)]).to(u.deg) ** 2
    hist, edges = np.histogram(theta_2, bins=100, range=u.Quantity([0, 30], u.deg**2))
    return hist, edges


def plot_theta2(theta_hist_dict):
    fig_theta2, ax_theta2 = plt.subplots(1, 1)
    alpha = 0.5
    cmap = plt.cm.get_cmap("viridis", len(theta_hist_dict))

    for i, (combiner, (hist, edges)) in enumerate(theta_hist_dict.items()):
        color = cmap(i)
        ax_theta2.stairs(
            hist,
            edges.to_value(u.deg**2),
            fill=False,
            alpha=alpha,
            color=color,
            label=combiner,
        )

    ax_theta2.set_title(r"$\theta^2$ for different StereoCombiner")
    ax_theta2.set_xlabel(r"$\theta^2$ / deg")
    ax_theta2.set_ylabel("# of events")
    ax_theta2.legend()
    ax_theta2.set_yscale("log")
    fig_theta2.tight_layout()

    return fig_theta2


def stack_lon_lat_hist(dl2_chunk):
    alt = dl2_chunk["disp_alt"]
    az = dl2_chunk["disp_az"]
    lon, lat = horizontal_to_telescope(
        alt,
        az,
        dl2_chunk["subarray_pointing_lat"],
        dl2_chunk["subarray_pointing_lon"],
    )

    lon_no_nans = lon[np.isfinite(lon)]
    lat_no_nans = lat[np.isfinite(lat)]
    # az_mod = (az_no_nans + 180) % 360 - 180

    hist, xedges, yedges = np.histogram2d(
        lon_no_nans.value, lat_no_nans.value, bins=(100, 100), range=[[-7, 7], [-7, 7]]
    )
    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]

    return hist, extent


def square_subplots(n, figsize_scale=4):
    cols = math.ceil(math.sqrt(n))
    rows = math.ceil(n / cols)
    fig, axes = plt.subplots(
        rows, cols, figsize=(cols * figsize_scale, rows * figsize_scale)
    )
    axes = np.atleast_1d(axes).flatten()
    return fig, axes


def plot_lon_lat(lon_lat_hist_dict, extent):
    n_combiner = len(lon_lat_hist_dict)
    fig, axes = square_subplots(n_combiner)

    for i, (combiner, hist) in enumerate(lon_lat_hist_dict.items()):
        ax = axes[i]
        h = ax.imshow(
            hist.T,
            origin="lower",
            extent=extent,
            cmap="viridis",
            norm=LogNorm(clip=True),
            aspect="equal",
        )
        ax.set_xlabel("Longitude / deg")
        ax.set_ylabel("Latitude / deg")
        ax.set_title(f"{combiner}")
        ax.set_aspect("equal")
        ax.set_xlim(-3, 3)
        ax.set_ylim(-3, 3)
        fig.colorbar(h, ax=ax)

    for ax in axes[n_combiner:]:
        ax.set_visible(False)

    fig.tight_layout()

    return fig


def plot_a_eff(irfs_dict):
    fig, ax = plt.subplots()
    for combiner, irfs_path in irfs_dict.items():
        aeff_table = EffectiveAreaTable2D.read(irfs_path, hdu="EFFECTIVE AREA")
        aeff_table.plot_energy_dependence(ax=ax, label=combiner, offset=[0.4] * u.deg)
    ax.legend()
    ax.set_yscale("log")
    ax.set_xlim(e_lim)
    ax.grid(which="both", linestyle=":")
    ax.set_title("Effective Area")
    ax.set_ylabel(r"Effective Area / $m^{2}$")
    ax.set_xlabel(r"$E_{True}$ / TeV")

    fig.tight_layout()
    return fig


def plot_angular_resolution(benchmarks_dict):
    fig, ax = plt.subplots()
    for combiner, benchmarks_path in benchmarks_dict.items():
        ang_res = QTable.read(benchmarks_path, hdu="ANGULAR RESOLUTION")
        ax.errorbar(
            (0.5 * (ang_res["ENERG_LO"] + ang_res["ENERG_HI"])).flatten(),
            ang_res["ANGULAR_RESOLUTION_68"].flatten(),
            xerr=0.5 * (ang_res["ENERG_HI"] - ang_res["ENERG_LO"]),
            ls="",
            label=f"{combiner}",
        )

    ax.set_xlim(e_lim)
    ax.set_ylim(0, 1)
    ax.set_title("Angular Resolution")
    ax.set_xscale("log")
    ax.set_xlabel(r"$E_{True}$ / TeV")
    ax.set_ylabel("68% Containment Angular Resolution / deg")
    ax.grid(which="both", linestyle=":")
    ax.legend()

    fig.tight_layout()
    return fig


def plot_energy_resolution(benchmarks_dict):
    fig, ax = plt.subplots()
    for combiner, benchmarks_path in benchmarks_dict.items():
        e_res = QTable.read(benchmarks_path, hdu="ENERGY BIAS RESOLUTION")
        ax.errorbar(
            (0.5 * (e_res["ENERG_LO"] + e_res["ENERG_HI"])).flatten(),
            e_res["RESOLUTION"].flatten(),
            xerr=0.5 * (e_res["ENERG_HI"] - e_res["ENERG_LO"]),
            ls="",
            label=f"{combiner}",
        )

    ax.set_xlim(e_lim)
    ax.set_ylim(0, 1)
    ax.set_title("Energy Resolution")
    ax.set_xscale("log")
    ax.set_xlabel(r"$E_{True}$ / TeV")
    ax.set_ylabel("Energy Resolution")
    ax.grid(which="both", linestyle=":")
    ax.legend()

    fig.tight_layout()
    return fig


def plot_sensitivity(benchmarks_dict):
    fig, ax = plt.subplots()
    for combiner, benchmarks_path in benchmarks_dict.items():
        sens = QTable.read(benchmarks_path, hdu="SENSITIVITY")
        ax.errorbar(
            (0.5 * (sens["ENERG_LO"] + sens["ENERG_HI"])).flatten(),
            sens["ENERGY_FLUX_SENSITIVITY"].flatten(),
            xerr=0.5 * (sens["ENERG_HI"] - sens["ENERG_LO"]),
            ls="",
            label=f"{combiner}",
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
    return fig
