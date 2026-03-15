from argparse import ArgumentParser
import matplotlib

from .irf_plots import (
    plot_angular_resolution,
    plots_cuts_distribution,
    plot_energy,
    plot_sensitivity,
    plot_a_eff,
)

if matplotlib.get_backend() == "pgf":
    from matplotlib.backends.backend_pgf import PdfPages
else:
    from matplotlib.backends.backend_pdf import PdfPages


parser = ArgumentParser()
parser.add_argument("--irfs-file", required=True)
parser.add_argument("--cuts-file", required=True)
parser.add_argument("--benchmark-file", required=True)
parser.add_argument("-o", "--output", required=True)
args = parser.parse_args()


def main(irfs_file, cuts_file, benchmark_file, output):
    fig_sens = plot_sensitivity(benchmark_file)
    fig_a_eff = plot_a_eff(irfs_file)
    fig_ang_res = plot_angular_resolution(benchmark_file)
    figs_energy = plot_energy(irfs_file, benchmark_file)
    fig_cuts = plots_cuts_distribution(cuts_file)
    # figs_bkg = plot_background_rate(irfs_file)

    with PdfPages(output) as pdf:
        pdf.savefig(fig_sens)
        pdf.savefig(fig_a_eff)
        pdf.savefig(fig_ang_res)
        for fig in figs_energy:
            pdf.savefig(fig)
        pdf.savefig(fig_cuts)


if __name__ == "__main__":
    main(**vars(args))
