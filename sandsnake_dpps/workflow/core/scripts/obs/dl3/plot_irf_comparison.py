import logging
from argparse import ArgumentParser
from contextlib import contextmanager

import astropy.units as u
import matplotlib
from gammapy.analysis import AnalysisConfig
from gammapy.data import DataStore
from matplotlib import pyplot as plt

from scriptutils.log import setup_logging

if matplotlib.get_backend() == "pgf":
    from matplotlib.backends.backend_pgf import PdfPages
else:
    from matplotlib.backends.backend_pdf import PdfPages

log = logging.getLogger(__name__)


# global units
@contextmanager
def autoscale_turned_off(ax=None):
    """
    https://stackoverflow.com/questions/38629830/how-to-turn-off-autoscaling-in-matplotlib-pyplot
    """
    ax = ax or plt.gca()
    lims = [ax.get_xlim(), ax.get_ylim()]
    yield
    ax.set_xlim(*lims[0])
    ax.set_ylim(*lims[1])


def mark_energy_region(ax, e_min, e_max):
    ax_min, ax_max = ax.get_xlim()
    if e_min is not None:
        ax.axvspan(ax_min, e_min, facecolor="black", alpha=0.3, hatch="/")
    if e_max is not None:
        ax.axvspan(ax_max, e_max, facecolor="black", alpha=0.3, hatch="/")


def title_irf_vs_offset(ax, energy):
    if matplotlib.rcParams["text.usetex"]:
        ax.set_title(
            r"$E \mathbin{=}  "
            + str(energy.to_value(u.GeV))
            + r"\, \unit{\giga\electronvolt}$",
        )
    else:
        ax.set_title(f"E = {energy}")


def title_irf_vs_energy(ax, offset):
    if matplotlib.rcParams["text.usetex"]:
        ax.set_title(
            r"$\text{Offset}:" + str(offset.to_value(u.deg)) + r" \unit{\degree}$",
        )
    else:
        ax.set_title(f"Offset: {offset}")


def x_vs_energy(ax):
    # log und label
    ax.set_xscale("log")
    if matplotlib.rcParams["text.usetex"]:
        ax.set_xlabel(r"$E_\text{reco} \mathbin{/} \unit{\giga\electronvolt}$")
    else:
        ax.set_xlabel("E_reco / GeV")


def x_vs_energy_true(ax):
    # log und label
    ax.set_xscale("log")
    if matplotlib.rcParams["text.usetex"]:
        ax.set_xlabel(r"$E_\text{true} \mathbin{/} \unit{\giga\electronvolt}$")
    else:
        ax.set_xlabel("E_true / GeV")


def x_vs_offset(ax):
    if matplotlib.rcParams["text.usetex"]:
        ax.set_xlabel(r"$\text{Offset} \mathbin{/} \unit{\degree}$")
    else:
        ax.set_xlabel("Offset / deg")


def y_bkg(ax):
    if matplotlib.rcParams["text.usetex"]:
        ax.set_ylabel(
            r"$\text{event rate} \mathord{/}"
            + r"\unit{\per\second\per\giga\electronvolt\per\steradian}$",
        )
    else:
        ax.set_ylabel("rate / s")


def y_edisp(ax):
    if matplotlib.rcParams["text.usetex"]:
        ax.set_ylabel(
            r"$edisp$",
        )
    else:
        ax.set_ylabel("edisp")


def y_psf(ax):
    if matplotlib.rcParams["text.usetex"]:
        ax.set_ylabel(
            r"$PSF \mathord{/} \unit{\degree}$",
        )
    else:
        ax.set_ylabel("PSF / deg")


def y_aeff(ax):
    ax.set_yscale("log")
    if matplotlib.rcParams["text.usetex"]:
        ax.set_ylabel(
            r"$A_\text{eff} \mathord{/} \unit{\meter\squared}$",
        )
    else:
        ax.set_ylabel("aeff / m^2")


def plot_psf_vs_energy(psf, ax, offset, **kwargs):
    energy_axis = psf.axes["energy_true"]
    e = energy_axis.center.to(u.GeV)
    radius = psf.containment_radius(
        energy_true=e,
        offset=offset,
        fraction=0.68,
    )

    ax.plot(e, radius, label=None, **kwargs)
    return ax


def plot_psf_vs_offset(psf, ax, energy, **kwargs):
    offset_axis = psf.axes["offset"]
    o = offset_axis.center
    radius = psf.containment_radius(
        energy_true=energy,
        offset=o,
        fraction=0.68,
    )

    ax.plot(o, radius, label=None, **kwargs)
    return ax


def plot_aeff_vs_energy(aeff, ax, offset, **kwargs):
    energy_axis = aeff.axes["energy_true"]
    e = energy_axis.center.to(u.GeV)
    area = aeff.evaluate(offset=offset, energy_true=e).to(u.m**2)
    ax.plot(e, area, label=None, **kwargs)
    return ax


def plot_aeff_vs_offset(aeff, ax, energy, **kwargs):
    offset_axis = aeff.axes["offset"]
    o = offset_axis.center
    area = aeff.evaluate(offset=o, energy_true=energy).to(u.m**2)
    ax.plot(o, area, label=None, **kwargs)
    return ax


def plot_edisp_vs_energy(edisp, ax, offset, **kwargs):
    energy_axis = edisp.axes["energy_true"]
    e = energy_axis.center.to(u.GeV)
    # add std? second axis?
    # -1 is invalid, all zeros
    kernel = edisp.to_edisp_kernel(offset=offset)
    bias = kernel.get_bias(e)
    ax.plot(e, bias.flatten(), label=None, **kwargs)
    return ax


def plot_edisp_vs_offset(edisp, ax, energy, **kwargs):
    offset_axis = edisp.axes["offset"]
    bias = []
    for o in offset_axis.center:
        kernel = edisp.to_edisp_kernel(offset=o)
        bias.append(kernel.get_bias(energy))
    ax.plot(offset_axis.center, bias, label=None, **kwargs)
    return ax


def plot_bkg2d_vs_energy(bkg2d, ax, offset, **kwargs):
    energy_axis = bkg2d.axes["energy"]
    e = energy_axis.center.to(u.GeV)
    rate = bkg2d.evaluate(offset=offset, energy=e).to(1 / (u.GeV * u.s * u.sr))
    ax.plot(e, rate, label=None, **kwargs)
    return ax


def plot_bkg2d_vs_offset(bkg2d, ax, energy, **kwargs):
    offset_axis = bkg2d.axes["offset"]
    o = offset_axis.center
    rate = bkg2d.evaluate(offset=o, energy=energy).to(1 / (u.GeV * u.s * u.sr))
    ax.plot(o, rate, label=None, **kwargs)
    return ax


def add_irf(  # noqa: PLR0913
    irf,
    axes,
    offset=0.4 * u.deg,
    energy=0.1 * u.TeV,
    setup=False,
    e_min=None,
    e_max=None,
):
    match irf.tag:
        case "edisp_2d":
            plot_edisp_vs_energy(irf, ax=axes[0], offset=offset, color="C0")
            plot_edisp_vs_offset(irf, ax=axes[1], energy=energy, color="C0")
            if setup:
                y_edisp(axes[0])
                x_vs_energy_true(axes[0])
        case "aeff_2d":
            plot_aeff_vs_energy(irf, ax=axes[0], offset=offset, color="C0")
            plot_aeff_vs_offset(irf, ax=axes[1], energy=energy, color="C0")
            if setup:
                y_aeff(axes[0])
                x_vs_energy_true(axes[0])
        case "psf_table":
            plot_psf_vs_energy(irf, ax=axes[0], offset=offset, color="C0")
            plot_psf_vs_offset(irf, ax=axes[1], energy=energy, color="C0")
            if setup:
                y_psf(axes[0])
                x_vs_energy_true(axes[0])
        case "bkg_2d":
            plot_bkg2d_vs_energy(irf, ax=axes[0], offset=offset, color="C0")
            plot_bkg2d_vs_offset(irf, ax=axes[1], energy=energy, color="C0")
            if setup:
                y_bkg(axes[0])
                x_vs_energy_true(axes[0])
        case "bkg_3d":
            plot_bkg2d_vs_energy(irf.to_2d(), ax=axes[0], offset=offset, color="C0")
            plot_bkg2d_vs_offset(irf.to_2d(), ax=axes[1], energy=energy, color="C0")
            if setup:
                y_bkg(axes[0])
                x_vs_energy_true(axes[0])
        case _:
            log.error(f"Could not handle irf: {irf.tag} ({type(irf)})")
    if setup:
        title_irf_vs_energy(axes[0], offset)
        title_irf_vs_offset(axes[1], energy)
        x_vs_offset(axes[1])
        axes[0].axvline(energy.to_value(u.GeV), ls="--")
        axes[1].axvline(offset.to_value(u.deg), ls="--")
        with autoscale_turned_off(ax=axes[0]):
            mark_energy_region(axes[0], e_min, e_max)


def main(input_path, config_path, output):
    ds = DataStore.from_file(input_path)

    config = AnalysisConfig.read(config_path)
    energy_axis_true_config = config.datasets.geom.axes.energy_true
    e_min = energy_axis_true_config.min.to_value(u.GeV)
    e_max = energy_axis_true_config.max.to_value(u.GeV)

    figs = {
        "psf": plt.subplots(1, 2, sharey=True),
        "aeff": plt.subplots(1, 2, sharey=True),
        "edisp": plt.subplots(1, 2, sharey=True),
        "bkg": plt.subplots(1, 2, sharey=True),
    }

    for j, o in enumerate(ds.get_observations()):
        for i in o.available_irfs:
            fig, axes = figs[i]
            add_irf(o.__getattribute__(i), axes, setup=j == 0, e_min=e_min, e_max=e_max)

    figures = []
    for fig, axes in figs.values():
        figures.append(fig)

    if output is None:
        plt.show()
    else:
        with PdfPages(output) as pdf:
            for fig in figures:
                pdf.savefig(fig)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-i", "--input-path", required=True)
    parser.add_argument("-o", "--output")
    parser.add_argument("--log-file")
    parser.add_argument("--config")  # the analysis config
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()
    setup_logging(logfile=args.log_file, verbose=args.verbose)

    main(args.input_path, args.config, args.output)
