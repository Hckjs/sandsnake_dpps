from argparse import ArgumentParser
import matplotlib

from core.scripts.mc.irf_plots import (
    plot_a_eff,
    plot_angular_resolution,
    plot_psf_radius,
    plots_cuts_distribution,
    plot_energy_resolution,
    plot_energy,
    plot_sensitivity,
    plot_background_rate_energy,
)

if matplotlib.get_backend() == "pgf":
    from matplotlib.backends.backend_pgf import PdfPages
else:
    from matplotlib.backends.backend_pdf import PdfPages


def main(irfs_file, cuts_file, benchmark_file, output):
    if irfs_file is None and cuts_file is None and benchmark_file is None:
        raise ValueError("No input files are given.")

    figs = []
    if benchmark_file is not None:
        figs.extend(plot_sensitivity(benchmark_file))
        figs.extend(plot_angular_resolution(benchmark_file))
        figs.extend(plot_energy_resolution(benchmark_file))

    if irfs_file is not None:
        figs.extend(plot_energy(irfs_file))
        figs.extend(plot_a_eff(irfs_file))
        figs.extend(plot_psf_radius(irfs_file))
        figs.extend(plot_background_rate_energy(irfs_file))

    if cuts_file is not None:
        figs.extend(plots_cuts_distribution(cuts_file))

    with PdfPages(output) as pdf:
        for fig in figs:
            pdf.savefig(fig)


def parse_args():
    parser = ArgumentParser()
    parser.add_argument("--irfs-file", required=False, default=None)
    parser.add_argument("--cuts-file", required=False, default=None)
    parser.add_argument("--benchmark-file", required=False, default=None)
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
