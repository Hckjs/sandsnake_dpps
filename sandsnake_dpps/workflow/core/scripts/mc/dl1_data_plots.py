from ctapipe.visualization import CameraDisplay
import numpy as np


## Histograms ##
def hist_total_intensity(data, ax, **kwargs):
    bins = np.logspace(0, 7, 30)
    ax.hist(
        data,
        bins=bins,
        stacked=False,
        histtype="step",
        density=True,
        fill=False,
        **kwargs,
    )
    ax.set_title("Hillas Intensity")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.legend()
    ax.set_xlabel("")


def hist_tels_with_trigger(data, ax, **kwargs):
    bins = 30
    ax.hist(
        data,
        bins=bins,
        stacked=False,
        histtype="step",
        density=True,
        fill=True,
        **kwargs,
    )
    ax.set_title("Triggered telescopes")
    ax.set_xlabel("")


def hist_trigger_per_tel(data, ax, **kwargs):
    bins = np.linspace(-0.5, 30.5, 31)
    ax.hist(
        data,
        bins=bins,
        stacked=False,
        histtype="step",
        density=True,
        fill=False,
        **kwargs,
    )
    ax.set_title("Tel type per subarray event")
    ax.set_xlabel("")
    ax.legend()


def hist_selected_pixels(data, ax, **kwargs):
    bins = np.logspace(0, 4, 30)
    ax.hist(
        data,
        bins=bins,
        stacked=False,
        histtype="step",
        density=True,
        fill=False,
        **kwargs,
    )
    ax.set_title("Selected pixels per tel event")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("")
    ax.legend()


def hist_identified_cherenkov_signal(data, ax, **kwargs):
    bins = 20
    ax.hist(
        data,
        bins=bins,
        stacked=False,
        histtype="step",
        density=True,
        fill=False,
        **kwargs,
    )
    ax.set_title("ICS per telescope event")
    ax.set_xlabel("")
    ax.legend()


def hist_tel_impact_distance(data, ax, **kwargs):
    bins = np.logspace(-1, 7, 30)
    ax.hist(
        data,
        bins=bins,
        stacked=False,
        histtype="step",
        density=True,
        fill=False,
        **kwargs,
    )
    ax.set_title("Telescope impact distance")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("")
    ax.legend()


# for diffuse
def hist_hillas_psi(data, ax, **kwargs):
    bins = 30
    ax.hist(
        data,
        bins=bins,
        stacked=False,
        histtype="step",
        density=True,
        fill=False,
        **kwargs,
    )
    ax.set_title("Hillas Psi")
    ax.set_xlabel("")
    ax.legend()


def hist_leakage(data, ax, **kwargs):
    bins = np.logspace(-1, 0, 30)
    ax.hist(
        data,
        bins=bins,
        stacked=False,
        histtype="step",
        density=True,
        fill=False,
        **kwargs,
    )
    ax.set_title("Leakage Intensity (width = 2)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("")
    ax.legend()


def hist_concentration_core(data, ax, **kwargs):
    bins = 30
    ax.hist(
        data,
        bins=bins,
        stacked=False,
        histtype="step",
        density=True,
        fill=False,
        **kwargs,
    )
    ax.set_title("Concentration Core")
    ax.set_xlabel("")
    ax.legend()


def hist_n_islands(data, ax, **kwargs):
    conc_data = np.concatenate(data).ravel()
    bins = np.linspace(-0.5, np.max(conc_data) + 0.5, np.max(conc_data) + 2)
    ax.hist(
        data,
        stacked=False,
        histtype="step",
        density=True,
        fill=False,
        bins=bins,
        **kwargs,
    )
    ax.set_title("Number of islands after cleanig")
    ax.set_yscale("log")
    ax.set_xlabel("")
    ax.legend()


## Energy dependent plots ##
def hillas_intensity_per_energy(dist1, dist2, fig, ax, tel, **kwargs):
    x_range = (np.nanmin(dist1), np.nanmax(dist1))
    y_range = (np.nanmin(dist2), np.nanmax(dist2))
    bins = [
        np.geomspace(x_range[0], x_range[1], 30),
        np.geomspace(y_range[0], y_range[1], 30),
    ]
    hist = ax.hist2d(
        dist1, dist2, bins=bins, norm="log", range=[x_range, y_range], **kwargs
    )
    fig.colorbar(hist[3], ax=ax)
    ax.set_title(f"{tel} Hillas Intensity")
    ax.set_xlabel("True Energy / TeV")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylabel("Hillas Intensity")


def ics_per_energy(dist1, dist2, fig, ax, tel, **kwargs):
    x_range = (np.nanmin(dist1), np.nanmax(dist1))
    y_range = (np.nanmin(dist2), np.nanmax(dist2))
    bins = [
        np.geomspace(x_range[0], x_range[1], 30),
        np.linspace(y_range[0], y_range[1], 30),
    ]
    hist = ax.hist2d(
        dist1, dist2, bins=bins, norm="log", range=[x_range, y_range], **kwargs
    )
    fig.colorbar(hist[3], ax=ax)
    ax.set_title(f"{tel} Identified Cherenkov Signal")
    ax.set_xlabel("True Energy / TeV")
    ax.set_xscale("log")
    ax.set_ylabel("ICS")


def hillas_wl_per_energy(dist1, dist2, fig, ax, tel, **kwargs):
    x_range = (np.nanmin(dist1), np.nanmax(dist1))
    y_range = (np.nanmin(dist2), np.nanmax(dist2))
    bins = [
        np.geomspace(x_range[0], x_range[1], 30),
        np.linspace(y_range[0], y_range[1], 30),
    ]
    hist = ax.hist2d(
        dist1, dist2, bins=bins, norm="log", range=[x_range, y_range], **kwargs
    )
    fig.colorbar(hist[3], ax=ax)
    ax.set_title(f"{tel} Hillas width/length")
    ax.set_xlabel("True Energy / TeV")
    ax.set_xscale("log")
    ax.set_ylabel("Hillas width/length")


def hillas_wl_true_per_energy(dist1, dist2, fig, ax, tel, **kwargs):
    x_range = (np.nanmin(dist1), np.nanmax(dist1))
    y_range = (np.nanmin(dist2), np.nanmax(dist2))
    bins = [
        np.geomspace(x_range[0], x_range[1], 30),
        np.linspace(y_range[0], y_range[1], 30),
    ]
    hist = ax.hist2d(
        dist1, dist2, bins=bins, norm="log", range=[x_range, y_range], **kwargs
    )
    fig.colorbar(hist[3], ax=ax)
    ax.set_title(f"{tel} Hillas true width/length")
    ax.set_xlabel("True Energy / TeV")
    ax.set_xscale("log")
    ax.set_ylabel("Hillas true width/length")


# for diffuse
def hillas_psi_per_energy(dist1, dist2, fig, ax, tel, **kwargs):
    x_range = (np.nanmin(dist1), np.nanmax(dist1))
    y_range = (np.nanmin(dist2), np.nanmax(dist2))
    bins = [
        np.geomspace(x_range[0], x_range[1], 30),
        np.linspace(y_range[0], y_range[1], 30),
    ]
    hist = ax.hist2d(
        dist1, dist2, bins=bins, norm="log", range=[x_range, y_range], **kwargs
    )
    fig.colorbar(hist[3], ax=ax)
    ax.set_title(f"{tel} Hillas Psi")
    ax.set_xlabel("True Energy / TeV")
    ax.set_xscale("log")
    ax.set_ylabel("Hillas Psi")


def intensity_max_per_energy(dist1, dist2, fig, ax, tel, **kwargs):
    x_range = (np.nanmin(dist1), np.nanmax(dist1))
    y_range = (np.nanmin(dist2), np.nanmax(dist2))
    bins = [
        np.geomspace(x_range[0], x_range[1], 30),
        np.geomspace(y_range[0], y_range[1], 30),
    ]
    hist = ax.hist2d(
        dist1, dist2, bins=bins, norm="log", range=[x_range, y_range], **kwargs
    )
    fig.colorbar(hist[3], ax=ax)
    ax.set_title(f"{tel} Intensity max")
    ax.set_xlabel("True Energy / TeV")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylabel("Intensity max")


# vergleich pro Node & particle type wäre noch interessant
def peak_time_max_min_per_energy(dist1, dist2, fig, ax, tel, **kwargs):
    x_range = (np.nanmin(dist1), np.nanmax(dist1))
    y_range = (np.nanmin(dist2), np.nanmax(dist2))
    bins = [
        np.geomspace(x_range[0], x_range[1], 30),
        np.linspace(0, 40, 30),
    ]
    hist = ax.hist2d(
        dist1, dist2, bins=bins, norm="log", range=[x_range, y_range], **kwargs
    )
    fig.colorbar(hist[3], ax=ax)
    ax.set_title(f"{tel} Peak time max-min")
    ax.set_xlabel("True Energy / TeV")
    ax.set_xscale("log")
    ax.set_ylabel("Peak time max-min")


def tel_impact_distance_per_energy(dist1, dist2, fig, ax, tel, **kwargs):
    x_range = (np.nanmin(dist1), np.nanmax(dist1))
    y_range = (np.nanmin(dist2), np.nanmax(dist2))
    bins = [
        np.geomspace(x_range[0], x_range[1], 30),
        np.geomspace(y_range[0], y_range[1], 30),
    ]
    hist = ax.hist2d(
        dist1, dist2, bins=bins, norm="log", range=[x_range, y_range], **kwargs
    )
    fig.colorbar(hist[3], ax=ax)
    ax.set_title(f"{tel} Telescope impact distance")
    ax.set_xlabel("True Energy / TeV")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylabel("Telescope impact distance")


## Reco vs. True ##
def hillas_intensity_true(dist1, dist2, fig, ax, tel, **kwargs):
    x_range = (np.nanmin(dist1), np.nanmax(dist1))
    y_range = (np.nanmin(dist2), np.nanmax(dist2))
    bins = [
        np.geomspace(x_range[0], x_range[1], 30),
        np.linspace(y_range[0], y_range[1], 30),
    ]
    hist = ax.hist2d(
        dist1, dist2, bins=bins, norm="log", range=[x_range, y_range], **kwargs
    )
    fig.colorbar(hist[3], ax=ax)
    ax.set_title(f"{tel} Hillas Intensity")
    ax.set_xlabel("True Hillas Intensity")
    ax.set_xscale("log")
    ax.set_ylabel("(Reco - True) / True")


def concentration_core_true(dist1, dist2, fig, ax, tel, **kwargs):
    x_range = (np.nanmin(dist1), np.nanmax(dist1))
    y_range = (np.nanmin(dist2), np.nanmax(dist2))
    bins = [
        np.linspace(x_range[0], x_range[1], 30),
        np.linspace(y_range[0], y_range[1], 30),
    ]
    hist = ax.hist2d(
        dist1, dist2, bins=bins, norm="log", range=[x_range, y_range], **kwargs
    )
    fig.colorbar(hist[3], ax=ax)
    ax.set_title(f"{tel} Concentration Core")
    ax.set_xlabel("True Concentration Core")
    ax.set_ylabel("(Reco - True) / True")


def n_islands_true(dist1, dist2, fig, ax, tel, **kwargs):
    x_range = (np.nanmin(dist1), np.nanmax(dist1))
    y_range = (np.nanmin(dist2), np.nanmax(dist2))
    bins = [
        np.linspace(x_range[0], x_range[1], 30),
        np.linspace(y_range[0], y_range[1], 30),
    ]
    hist = ax.hist2d(
        dist1, dist2, bins=bins, norm="log", range=[x_range, y_range], **kwargs
    )
    fig.colorbar(hist[3], ax=ax)
    ax.set_title(f"{tel} n Islands")
    ax.set_xlabel("True n Islands")
    ax.set_ylabel("(Reco - True) / True")


def tel_impact_distance_true(dist1, dist2, fig, ax, tel, **kwargs):
    x_range = (np.nanmin(dist1), np.nanmax(dist1))
    y_range = (np.nanmin(dist2), np.nanmax(dist2))
    bins = [
        np.linspace(0, 1200, 30),
        np.linspace(0, 1200, 30),
    ]
    hist = ax.hist2d(
        dist1, dist2, bins=bins, norm="log", range=[x_range, y_range], **kwargs
    )
    fig.colorbar(hist[3], ax=ax)
    ax.set_title(f"{tel} Telescope impact distance")
    ax.set_xlabel("True Telescope impact distance")
    ax.set_ylabel("(Reco - True) / True")


## Camera Display ##
def cam_intensity_mean(geom, image, ax, **kwargs):
    disp = CameraDisplay(geom, image=image, ax=ax, **kwargs)
    disp.cmap = "viridis"
    disp.add_colorbar(ax=ax)
    ax.set_title(f"Mean {geom.name}")


def cam_intensity_std(geom, image, ax, **kwargs):
    disp = CameraDisplay(geom, image=image, ax=ax, **kwargs)
    disp.cmap = "viridis"
    disp.add_colorbar(ax=ax)
    ax.set_title(f"Std {geom.name}")
