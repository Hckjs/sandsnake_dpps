import astropy.units as u
from astropy.table import QTable
from astropy.io import fits
import matplotlib.pyplot as plt
from pyirf.spectral import CRAB_MAGIC_JHEAP2015
import numpy as np
from pathlib import Path

from ctapipe.irf.spectra import Spectra, SPECTRA, ENERGY_FLUX_UNIT
from gammapy.irf import EnergyDispersion2D, EffectiveAreaTable2D, Background2D, PSF3D
from common.plotting.colors import CTAO_COLORS, CTAO_CMAP_R

REPO_ROOT = Path(__file__).resolve().parents[5]


def plots_cuts_distribution(cuts_path):
    figs = []
    with fits.open(cuts_path) as hdul:
        hdu_names = [hdu.name for hdu in hdul]

    # G/H
    gh_cut = QTable.read(cuts_path, hdu="GH_CUTS")
    fig_gh, ax_gh = plt.subplots()

    e_lim = [gh_cut["low"][0].to_value(u.TeV), gh_cut["high"][-1].to_value(u.TeV)]

    ax_gh.errorbar(
        gh_cut["center"],
        gh_cut["cut"],
        xerr=(gh_cut["center"] - gh_cut["low"], gh_cut["high"] - gh_cut["center"]),
        ls="",
    )

    ax_gh.set_xlim(e_lim)
    ax_gh.set_xscale("log")
    ax_gh.set_title(r"$\gamma$/H Cuts")
    ax_gh.set_ylabel(r"$\gamma$/H Cut")
    ax_gh.set_xlabel(r"$E_{Reco}$ / TeV")
    ax_gh.grid(which="both")
    figs.append(fig_gh)

    # Theta
    if "RAD_MAX" in hdu_names:
        theta_cut = QTable.read(cuts_path, hdu="RAD_MAX")
        fig_theta, ax_theta = plt.subplots()

        ax_theta.errorbar(
            theta_cut["center"],
            theta_cut["cut"],
            xerr=(
                theta_cut["center"] - theta_cut["low"],
                theta_cut["high"] - theta_cut["center"],
            ),
            ls="",
        )
        ax_theta.set_xlim(e_lim)
        ax_theta.set_xscale("log")
        ax_theta.set_title(r"$\Theta$ Cuts")
        ax_theta.set_ylabel(r"$\Theta$ Cut")
        ax_theta.set_xlabel(r"$E_{Reco}$ / TeV")
        ax_theta.grid(which="both")
        figs.append(fig_theta)

    # Multiplicity
    if "MULTIPLICITY_CUTS" in hdu_names:
        mult_cut = QTable.read(cuts_path, hdu="MULTIPLICITY_CUTS")
        fig_mult, ax_mult = plt.subplots()

        ax_mult.errorbar(
            mult_cut["center"],
            mult_cut["cut"],
            xerr=(
                mult_cut["center"] - mult_cut["low"],
                mult_cut["high"] - mult_cut["center"],
            ),
            ls="",
        )
        ax_mult.set_xlim(e_lim)
        ax_mult.set_xscale("log")
        ax_mult.set_title(r"Event Multiplicity Cuts")
        ax_mult.set_ylabel(r"$N_\mathrm{Multiplicity}$ Cut")
        ax_mult.set_xlabel(r"$E_{Reco}$ / TeV")
        ax_mult.grid(which="both")
        figs.append(fig_mult)

    return figs


def plot_a_eff(irfs_path):
    fig, ax = plt.subplots()
    aeff_table = EffectiveAreaTable2D.read(irfs_path, hdu="EFFECTIVE AREA")

    aeff_table.plot_energy_dependence(
        ax=ax, offset=aeff_table.axes["offset"].center, alpha=0.7
    )

    ax.set_yscale("log")
    ax.set_xlim(
        aeff_table.axes["energy_true"].edges[0],
        aeff_table.axes["energy_true"].edges[-1],
    )
    ax.grid(which="both")
    ax.set_title("Effective Area")
    ax.set_ylabel("Effective Area / $m^{2}$")
    ax.set_xlabel(r"$E_{True}$ / TeV")
    return fig


def plot_energy(irfs_path, benchmarks_path):
    figs = []
    edisp_table = EnergyDispersion2D.read(irfs_path, hdu="ENERGY DISPERSION")
    e_resolution = QTable.read(benchmarks_path, hdu="ENERGY BIAS RESOLUTION")

    fig_res, ax_res = plt.subplots(1, 1)
    fig_bias, ax_bias = plt.subplots(1, 1)
    fig_mat, ax_mat = plt.subplots(1, 1)

    # Energy Resolution
    fov_centers = (
        0.5 * (e_resolution["THETA_LO"] + e_resolution["THETA_HI"])
    ).flatten()
    for i, offset in enumerate(fov_centers):
        ax_res.errorbar(
            (0.5 * (e_resolution["ENERG_LO"] + e_resolution["ENERG_HI"])).flatten(),
            e_resolution["RESOLUTION"][:, i, :].flatten(),
            xerr=0.5 * (e_resolution["ENERG_HI"] - e_resolution["ENERG_LO"]),
            ls="",
            alpha=0.7,
            label=f"offset = {offset}",
        )

    ax_res.set_xscale("log")
    ax_res.set_xlim(
        e_resolution["ENERG_LO"][0][0].to_value(u.TeV),
        e_resolution["ENERG_HI"][0][-1].to_value(u.TeV),
    )
    ax_res.set_title("Energy Resolution")
    ax_res.set_xlabel(r"$E_{True}$ / TeV")
    ax_res.set_ylabel("Energy Resolution")
    ax_res.grid(which="both")
    ax_res.legend()

    # Energy Bias
    edisp_table.plot_bias(
        ax=ax_bias, offset=fov_centers[0], add_cbar=True, cmap=CTAO_CMAP_R
    )
    ax_bias.set_title(f"Energy Bias at Offset = {fov_centers[0]:.1f}")
    ax_bias.set_xlabel(r"$E_{True}$ / TeV")
    ax_bias.set_ylabel(r"$E_{Reco}/E_{True}$")

    # Energy Migration Matrix
    x = np.linspace(0, 1000, 10)
    edisp_kernel = edisp_table.to_edisp_kernel(offset=fov_centers[0])
    edisp_kernel.plot_matrix(ax=ax_mat, add_cbar=True, cmap=CTAO_CMAP_R)

    ax_mat.plot(x, x, color="black", linestyle="--")
    ax_mat.set_title(f"Energy Migration Matrix at Offset = {fov_centers[0]:.1f}")
    ax_mat.set_xlabel(r"$E_{True}$ / TeV")
    ax_mat.set_ylabel(r"$E_{Reco}$ / TeV")

    figs.append(fig_res)
    figs.append(fig_mat)
    figs.append(fig_bias)
    return figs


def plot_background_rate_energy(irfs_path):
    """Plot background rate as function of reconstructed energy for one offset."""
    bkg_table = Background2D.read(irfs_path, hdu="BACKGROUND")
    fig, ax = plt.subplots()

    bkg_table.plot_energy_dependence(
        ax=ax, offset=bkg_table.axes["offset"].center, alpha=0.7
    )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_title("Background Rate")
    ax.set_xlabel(r"$E_\mathrm{Reco}$ / TeV")
    ax.set_ylabel(f"Background rate / {bkg_table.unit}")
    ax.grid(which="both")
    return fig


def plot_psf_radius(irfs_path):
    psf = PSF3D.read(irfs_path, hdu="PSF")
    fig, ax = plt.subplots()

    psf.plot_containment_radius(ax=ax, add_cbar=True, cmap=CTAO_CMAP_R)

    ax.set_xscale("log")
    ax.set_title("PSF 68% Containment Radius")
    ax.set_xlabel(r"$E_\mathrm{Reco}$ / TeV")
    ax.set_ylabel(r"FoV Offset / $^\circ$")
    ax.grid(which="both")
    return fig


def plot_angular_resolution(benchmarks_path):
    ang_res = QTable.read(benchmarks_path, hdu="ANGULAR RESOLUTION")
    # Compare with requirements
    fig_comp, ax_comp = plt.subplots()

    fov_centers = (0.5 * (ang_res["THETA_LO"] + ang_res["THETA_HI"])).flatten()
    e_type = ang_res["E_TYPE"]
    e_lim = [
        ang_res["ENERG_LO"][0][0].to_value(u.TeV),
        ang_res["ENERG_HI"][0][-1].to_value(u.TeV),
    ]

    ax_comp.errorbar(
        (0.5 * (ang_res["ENERG_LO"] + ang_res["ENERG_HI"])).flatten(),
        ang_res["ANGULAR_RESOLUTION_68"][:, 0, :].flatten(),
        xerr=0.5 * (ang_res["ENERG_HI"] - ang_res["ENERG_LO"]),
        ls="",
        label="R68",
    )

    resources_path = REPO_ROOT / "resources"
    ctao_req_e, ctao_req_ang = np.loadtxt(
        resources_path / "CTA_Requirements/cta_requirements_North-50h-AngRes.dat",
        unpack=True,
    )
    ax_comp.plot(
        ctao_req_e,
        ctao_req_ang,
        color=CTAO_COLORS["interstellar_indigo"],
        label="CTAO-N Req. (50h)",
        alpha=0.8,
    )

    ax_comp.set_xlim(e_lim)
    ax_comp.set_ylim(0, 0.3)
    ax_comp.set_xscale("log")
    ax_comp.set_title(f"Angular Resolution at Offset = {fov_centers[0]:.1f}")
    ax_comp.set_xlabel(rf"$E_{{{e_type}}}$ / TeV")
    ax_comp.set_ylabel(r"Angular Resolution / $^{\circ}$")
    ax_comp.grid(which="both")
    ax_comp.legend()

    # Plot R68 for all FoV bins
    fig_fov, ax_fov = plt.subplots()
    for i, offset in enumerate(fov_centers):
        ax_fov.errorbar(
            (0.5 * (ang_res["ENERG_LO"] + ang_res["ENERG_HI"])).flatten(),
            ang_res["ANGULAR_RESOLUTION_68"][:, i, :].flatten(),
            xerr=0.5 * (ang_res["ENERG_HI"] - ang_res["ENERG_LO"]),
            ls="",
            alpha=0.7,
            label=f"offset = {offset}",
        )

    ax_fov.set_xlim(e_lim)
    ax_fov.set_xscale("log")
    ax_fov.set_title("Angular Resolution")
    ax_fov.set_xlabel(rf"$E_{{{e_type}}}$ / TeV")
    ax_fov.set_ylabel(r"68% Angular Resolution / $^{\circ}$")
    ax_fov.grid(which="both")
    ax_fov.legend()

    # Plot all quantiles present in file
    fig, ax = plt.subplots()

    for col in [c for c in ang_res.colnames if c.startswith("ANGULAR_RESOLUTION")]:
        ax.errorbar(
            (0.5 * (ang_res["ENERG_LO"] + ang_res["ENERG_HI"])).flatten(),
            ang_res[col][:, 0, :].flatten(),
            xerr=0.5 * (ang_res["ENERG_HI"] - ang_res["ENERG_LO"]),
            ls="",
            label=f"R{col[-2:]}",
        )

    ax.set_xlim(e_lim)
    ax.set_xscale("log")
    ax.set_title(f"Angular Resolution at Offset = {fov_centers[0]:.1f}")
    ax.set_xlabel(rf"$E_{{{e_type}}}$ / TeV")
    ax.set_ylabel(r"Angular Resolution / $^{\circ}$")
    ax.grid(which="both")
    ax.legend()
    return [fig_comp, fig_fov, fig]


def plot_sensitivity(benchmarks_path):
    sens = QTable.read(benchmarks_path, hdu="SENSITIVITY")
    fig, ax = plt.subplots()

    fov_centers = (0.5 * (sens["THETA_LO"] + sens["THETA_HI"])).flatten()
    e_lim = [
        sens["ENERG_LO"][0][0].to_value(u.TeV),
        sens["ENERG_HI"][0][-1].to_value(u.TeV),
    ]

    ax.errorbar(
        (0.5 * (sens["ENERG_LO"] + sens["ENERG_HI"])).flatten(),
        sens["ENERGY_FLUX_SENSITIVITY"][:, 0, :].flatten(),
        xerr=0.5 * (sens["ENERG_HI"] - sens["ENERG_LO"]),
        ls="",
    )

    add_sensitivity_comparisons(ax, energy_limits=e_lim)

    ax.set_ylim(3.0e-14, 1.0e-9)
    ax.set_xlim(e_lim)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_title(f"Sensitivity at Offset = {fov_centers[0]:.1f}")
    ax.set_xlabel(r"$E_{Reco}$ / TeV")
    ax.set_ylabel(r"$E^{2} \times$ Flux Sensitivity / $erg \cdot cm^{-2} \cdot s^{-1}$")
    ax.grid(which="both")
    ax.legend(loc="upper right", fontsize="small")

    # Plot sensitivity for all FoV bins
    fig_fov, ax_fov = plt.subplots()

    for i, offset in enumerate(fov_centers):
        ax_fov.errorbar(
            (0.5 * (sens["ENERG_LO"] + sens["ENERG_HI"])).flatten(),
            sens["ENERGY_FLUX_SENSITIVITY"][:, i, :].flatten(),
            xerr=0.5 * (sens["ENERG_HI"] - sens["ENERG_LO"]),
            ls="",
            alpha=0.7,
            label=f"offset = {offset}",
        )

    ax_fov.set_xlim(e_lim)
    ax_fov.set_xscale("log")
    ax_fov.set_yscale("log")
    ax_fov.set_title("Sensitivity")
    ax_fov.set_xlabel(r"$E_{Reco}$ / TeV")
    ax_fov.set_ylabel(
        r"$E^{2} \times$ Flux Sensitivity / $erg \cdot cm^{-2} \cdot s^{-1}$"
    )
    ax_fov.grid(which="both")
    ax_fov.legend()
    return [fig, fig_fov]


def add_sensitivity_comparisons(ax, energy_limits):
    # Plot Crab SED
    plot_Crab_SED(
        energy_limits[0] * u.TeV,
        energy_limits[1] * u.TeV,
        percentage=100,
        ax=ax,
        label="100% Crab",
        alpha=0.5,
    )
    plot_Crab_SED(
        energy_limits[0] * u.TeV,
        energy_limits[1] * u.TeV,
        percentage=10,
        ax=ax,
        linestyle="--",
        label="10% Crab",
        alpha=0.5,
    )
    plot_Crab_SED(
        energy_limits[0] * u.TeV,
        energy_limits[1] * u.TeV,
        percentage=1,
        ax=ax,
        linestyle=":",
        label="1% Crab",
        alpha=0.5,
    )

    resources_path = REPO_ROOT / "resources"
    ctao_req_e, ctao_req_sens = np.loadtxt(
        resources_path / "CTA_Requirements/cta_requirements_North-50h.dat",
        unpack=True,
    )
    ax.plot(
        ctao_req_e,
        ctao_req_sens,
        color=CTAO_COLORS["interstellar_indigo"],
        label="CTAO-N Req. (50h)",
        alpha=0.8,
    )

    veritas_data = np.loadtxt(
        resources_path / "VERITAS/VERITAS_V6_std_50hr_5sigma_VERITAS2014_DiffSens.dat",
    )
    veritas_e = veritas_data[:, 0]
    veritas_rel_flux = veritas_data[:, 1]
    ax.plot(
        veritas_e,
        (
            veritas_rel_flux
            * SPECTRA[Spectra.CRAB_HEGRA](veritas_e * u.TeV)
            * (veritas_e * u.TeV) ** 2
        ).to(ENERGY_FLUX_UNIT),
        color=CTAO_COLORS["cosmic_azure"],
        label="VERITAS (50h)",
        linestyle="dashed",
        alpha=0.8,
    )

    magic_e, magic_flux = np.genfromtxt(
        resources_path / "MAGIC/MAGIC_differential_sensitivity_50hr_2025.csv",
        unpack=True,
        skip_header=2,
    )
    ax.plot(
        (magic_e * u.GeV).to(u.TeV).value,
        magic_flux * u.Unit("TeV cm-2 s-1").to(ENERGY_FLUX_UNIT),
        color=CTAO_COLORS["cosmic_azure"],
        label="MAGIC (50h)",
        linestyle="dotted",
        alpha=0.8,
    )


def plot_Crab_SED(emin, emax, percentage=100, ax=None, **kwargs):
    """
    Plot a percentage of the Crab SED

    Parameters
    ----------
    emin: `astropy.units.quantity.Quantity` compatible with energies
    emax:  astropy.units.quantity.Quantity compatible with energies
    percentage:  `float`  percentage of the Crab Nebula to be plotted
    ax:    `matplotlib.pyplot.axis`
    kwargs: kwargs for `matplotlib.pyplot.plot`

    Returns
    -------
    ax:    `matplotlib.pyplot.axis`
    """
    ax = plt.gca() if ax is None else ax

    energy = np.geomspace(emin.to(u.TeV), emax.to(u.TeV), 40)

    if percentage == 100:
        kwargs.setdefault("label", "Crab (MAGIC JHEAP 2015)")
    else:
        kwargs.setdefault("label", f"{percentage}% Crab (MAGIC JHEAP 2015)")

    kwargs.setdefault("color", "gray")

    ax.plot(
        energy.value,
        (percentage / 100.0 * (energy**2 * CRAB_MAGIC_JHEAP2015(energy)))
        .to(u.erg / (u.s * u.cm**2))
        .value,
        **kwargs,
    )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Energy / TeV")
    ax.set_ylabel(r"E$^2$ $\frac{\mathrm{dN}}{\mathrm{dE}}$ / erg cm$^{-2}$ s$^{-1}$")
    ax.legend()
    return ax
