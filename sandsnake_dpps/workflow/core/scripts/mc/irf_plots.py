import astropy.units as u
from astropy.table import QTable
import matplotlib.pyplot as plt
from pyirf.spectral import CRAB_MAGIC_JHEAP2015
import numpy as np
from pathlib import Path

from ctapipe.irf.spectra import Spectra, SPECTRA, ENERGY_FLUX_UNIT
from gammapy.irf import EnergyDispersion2D, EffectiveAreaTable2D, Background2D

REPO_ROOT = Path(__file__).resolve().parents[5]

offset = 0.7 * u.deg
e_lim = [5.0e-3, 5.0e2]


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


def plot_a_eff(irfs_path):
    fig, ax = plt.subplots()
    aeff_table = EffectiveAreaTable2D.read(irfs_path, hdu="EFFECTIVE AREA")
    aeff_table.plot_energy_dependence(ax=ax, offset=[offset])
    ax.set_yscale("log")
    ax.set_xlim(e_lim)
    ax.grid(which="both", linestyle=":")
    ax.set_title("Effective Area")
    ax.set_ylabel("Effective Area / $m^{2}$")
    ax.set_xlabel(r"$E_{True}$ / TeV")
    ax.legend().remove()

    fig.tight_layout()
    return fig


def plot_energy(irfs_path, benchmarks_path):
    figs = []
    edisp_table = EnergyDispersion2D.read(irfs_path, hdu="ENERGY DISPERSION")
    e_resolution = QTable.read(benchmarks_path, hdu="ENERGY BIAS RESOLUTION")

    fig_res, ax_res = plt.subplots(1, 1)
    fig_bias, ax_bias = plt.subplots(1, 1)
    fig_mat, ax_mat = plt.subplots(1, 1)

    # Energy Resolution
    ax_res.errorbar(
        (0.5 * (e_resolution["ENERG_LO"] + e_resolution["ENERG_HI"])).flatten(),
        e_resolution["RESOLUTION"].flatten(),
        xerr=0.5 * (e_resolution["ENERG_HI"] - e_resolution["ENERG_LO"]),
        ls="",
    )

    ax_res.set_xscale("log")
    ax_res.set_xlim(e_lim)
    ax_res.set_title("Energy Resolution")
    ax_res.set_xlabel(r"$E_{True}$ / TeV")
    ax_res.set_ylabel("Energy Resolution")
    ax_res.grid(which="both", linestyle=":")
    fig_res.tight_layout()

    # Energy Bias
    edisp_table.plot_bias(ax=ax_bias, offset=offset, add_cbar=True)
    ax_bias.set_title("Energy Bias")
    ax_bias.set_xlabel(r"$E_{True}$ / TeV")
    ax_bias.set_ylabel(r"$E_{Reco}/E_{True}$")
    fig_bias.tight_layout()

    # Energy Migration Matrix
    x = np.linspace(0, 1000, 10)
    edisp_kernel = edisp_table.to_edisp_kernel(offset=offset)
    edisp_kernel.plot_matrix(
        ax=ax_mat,
        add_cbar=True,
    )
    ax_mat.plot(x, x, color="black", linestyle="--")
    ax_mat.set_title("Energy Migration Matrix")
    ax_mat.set_xlabel(r"$E_{True}$ / TeV")
    ax_mat.set_ylabel(r"$E_{Reco}$ / TeV")
    fig_mat.tight_layout()

    figs.append(fig_res)
    figs.append(fig_mat)
    figs.append(fig_bias)

    return figs


def plot_angular_resolution(benchmarks_path):
    # TODO:25/50/95 auch mit rein?
    ang_res = QTable.read(benchmarks_path, hdu="ANGULAR RESOLUTION")
    fig, ax = plt.subplots()

    ax.errorbar(
        (0.5 * (ang_res["ENERG_LO"] + ang_res["ENERG_HI"])).flatten(),
        ang_res["ANGULAR_RESOLUTION_68"].flatten(),
        xerr=0.5 * (ang_res["ENERG_HI"] - ang_res["ENERG_LO"]),
        ls="",
    )

    ax.set_xscale("log")
    ax.set_xlim(e_lim)
    ax.set_title("Angular Resolution")
    ax.set_xlabel(r"$E_{True}$ / TeV")
    ax.set_ylabel("68% Containment Angular Resolution / deg")
    ax.grid(which="both", linestyle=":")

    fig.tight_layout()
    return fig


def plots_cuts_distribution(cuts_path):
    # cuts distribution
    theta_cut = QTable.read(cuts_path, hdu="RAD_MAX")
    gh_cut = QTable.read(cuts_path, hdu="GH_CUTS")

    fig_cuts, ax_cuts = plt.subplots(2, 1)

    # Theta
    ax_cuts[0].errorbar(
        theta_cut["center"],
        theta_cut["cut"],
        xerr=(
            theta_cut["center"] - theta_cut["low"],
            theta_cut["high"] - theta_cut["center"],
        ),
    )
    ax_cuts[0].set_xlim(e_lim)
    ax_cuts[0].set_xscale("log")
    ax_cuts[0].set_title(r"$\Theta$ cuts")
    ax_cuts[0].set_ylabel(r"$\Theta$ cuts")
    ax_cuts[0].set_xlabel(r"$E_{True}$ / TeV")
    ax_cuts[0].grid(which="both", linestyle=":")

    # G/H
    ax_cuts[1].errorbar(
        gh_cut["center"],
        gh_cut["cut"],
        xerr=(gh_cut["center"] - gh_cut["low"], gh_cut["high"] - gh_cut["center"]),
    )

    ax_cuts[1].set_xlim(e_lim)
    ax_cuts[1].set_xscale("log")
    ax_cuts[1].set_title(r"$\gamma$/H Cuts")
    ax_cuts[1].set_ylabel(r"$\gamma$/H Cuts")
    ax_cuts[1].set_xlabel(r"$E_{True}$ / TeV")
    ax_cuts[1].grid(which="both", linestyle=":")

    fig_cuts.tight_layout()
    return fig_cuts


def plot_background_rate(irfs_path):
    figs = []
    # background rate
    bkg_table = Background2D.read(irfs_path, hdu="BACKGROUND")
    fig_off_dep, axes_off_dep = plt.subplots(1, 1)
    fig_bkg, axes_bkg = plt.subplots(1, 1)
    fig_e_dep, axes_e_dep = plt.subplots(1, 1)
    fig_spec, axes_spec = plt.subplots(1, 1)
    true_energy = [0.1, 0.5, 1, 10] * u.TeV

    bkg_table.plot_offset_dependence(ax=axes_off_dep, energy=true_energy)
    axes_off_dep.grid(which="both")
    axes_off_dep.set_title("Background rate - offset dependence")
    fig_off_dep.tight_layout()
    figs.append(fig_off_dep)

    bkg_table.plot(ax=axes_bkg)
    axes_bkg.grid()
    axes_bkg.set_title("Background rate - energy offset dependence")
    fig_bkg.tight_layout()
    figs.append(fig_bkg)

    bkg_offset = [0, 1, 2, 3, 4, 5] * u.deg
    labels = []
    for o in bkg_offset:
        bkg_table.plot_energy_dependence(ax=axes_e_dep, offset=[o])
        labels.append(f"offset = {o}")
    axes_e_dep.legend().remove()
    axes_e_dep.legend(labels)
    axes_e_dep.grid(which="both")
    axes_e_dep.set_title("For different offsets")
    fig_e_dep.tight_layout()
    figs.append(fig_e_dep)

    bkg_table.plot_spectrum(ax=axes_spec)
    axes_spec.grid(which="both")
    axes_spec.legend().remove()
    axes_spec.set_title("Integrated spectrum")
    fig_spec.tight_layout()
    figs.append(fig_spec)

    return figs


def plot_sensitivity(benchmarks_path):
    sens = QTable.read(benchmarks_path, hdu="SENSITIVITY")
    fig, ax = plt.subplots()

    ax.errorbar(
        (0.5 * (sens["ENERG_LO"] + sens["ENERG_HI"])).flatten(),
        sens["ENERGY_FLUX_SENSITIVITY"].flatten(),
        xerr=0.5 * (sens["ENERG_HI"] - sens["ENERG_LO"]),
        ls="",
        label=f"Offset: {offset:.1f}",
    )

    add_sensitivity_comparisons(ax)

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


def add_sensitivity_comparisons(ax, energy_limits=e_lim):
    # Plot Crab SED
    plot_Crab_SED(
        energy_limits[0] * u.TeV,
        energy_limits[1] * u.TeV,
        percentage=100,
        ax=ax,
        label="100% Crab",
    )
    plot_Crab_SED(
        energy_limits[0] * u.TeV,
        energy_limits[1] * u.TeV,
        percentage=10,
        ax=ax,
        linestyle="--",
        label="10% Crab",
    )
    plot_Crab_SED(
        energy_limits[0] * u.TeV,
        energy_limits[1] * u.TeV,
        percentage=1,
        ax=ax,
        linestyle=":",
        label="1% Crab",
    )

    resources_path = REPO_ROOT / "resources"
    ctao_req_e, ctao_req_sens = np.loadtxt(
        resources_path / "CTA_Requirements/cta_requirements_North-50h.dat",
        unpack=True,
    )
    ax.plot(
        ctao_req_e,
        ctao_req_sens,
        color="red",
        label="CTAO North Requirement (50h)",
        alpha=0.5,
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
        color="blue",
        label="VERITAS (50h)",
        linestyle="dashed",
        alpha=0.5,
    )

    magic_e, magic_flux = np.genfromtxt(
        resources_path / "MAGIC/MAGIC_differential_sensitivity_50hr_2025.csv",
        unpack=True,
        skip_header=2,
    )
    ax.plot(
        (magic_e * u.GeV).to(u.TeV).value,
        magic_flux * u.Unit("TeV cm-2 s-1").to(ENERGY_FLUX_UNIT),
        color="green",
        label="MAGIC (50h)",
        linestyle="dashed",
        alpha=0.5,
    )
