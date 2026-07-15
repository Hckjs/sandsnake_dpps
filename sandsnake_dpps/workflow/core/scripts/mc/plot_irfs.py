from argparse import ArgumentParser
import matplotlib

from core.scripts.mc.irf_plots import (
    plot_angular_resolution,
    plots_cuts_distribution,
    plot_energy,
    plot_sensitivity,
    plot_a_eff,
    plot_background_rate_energy,
    plot_psf_radius,
)

if matplotlib.get_backend() == "pgf":
    from matplotlib.backends.backend_pgf import PdfPages
else:
    from matplotlib.backends.backend_pdf import PdfPages


def main(irfs_file, cuts_file, benchmark_file, output):
    fig_sens = plot_sensitivity(benchmark_file)
    fig_a_eff = plot_a_eff(irfs_file)
    figs_ang_res = plot_angular_resolution(benchmark_file)
    figs_energy = plot_energy(irfs_file, benchmark_file)
    figs_cuts = plots_cuts_distribution(cuts_file)
    fig_bkg = plot_background_rate_energy(irfs_file)
    fig_psf = plot_psf_radius(irfs_file)

    with PdfPages(output) as pdf:
        pdf.savefig(fig_sens)
        for fig in figs_ang_res:
            pdf.savefig(fig)
        pdf.savefig(fig_a_eff)
        pdf.savefig(fig_psf)
        for fig in figs_energy + figs_cuts:
            pdf.savefig(fig)
        pdf.savefig(fig_bkg)


def parse_args():
    parser = ArgumentParser()
    parser.add_argument("--irfs-file", required=True)
    parser.add_argument("--cuts-file", required=True)
    parser.add_argument("--benchmark-file", required=True)
    parser.add_argument("-o", "--output", required=True)
    return parser.parse_args()


def main_from_snakemake(snakemake):
    main(
        snakemake.input.irfs,
        snakemake.input.cuts,
        snakemake.input.benchmarks,
        snakemake.output[0],
    )


if "snakemake" in globals():
    main_from_snakemake(snakemake)  # noqa: F821
elif __name__ == "__main__":
    args = parse_args()
    main(args.irfs_file, args.cuts_file, args.benchmark_file, args.output)
