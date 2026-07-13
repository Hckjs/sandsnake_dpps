from argparse import ArgumentParser
from pathlib import Path

from plugins.stereo_combiner.scripts.stereo_combiner_plots import (
    plot_a_eff,
    plot_angular_resolution,
    plot_sensitivity,
    plot_energy_resolution,
)

import matplotlib

if matplotlib.get_backend() == "pgf":
    from matplotlib.backends.backend_pgf import PdfPages
else:
    from matplotlib.backends.backend_pdf import PdfPages


def create_file_dict(file_list):
    file_dict = {}
    for file in file_list:
        combiner = Path(file).stem.split("_")[-1].split(".")[0]
        file_dict[combiner] = file

    return file_dict


def main(irf_files, benchmark_files, output):
    irf_file_dict = create_file_dict(irf_files)
    benchmark_file_dict = create_file_dict(benchmark_files)

    fig_sens = plot_sensitivity(benchmark_file_dict)
    fig_a_eff = plot_a_eff(irf_file_dict)
    fig_ang_res = plot_angular_resolution(benchmark_file_dict)
    fig_e_res = plot_energy_resolution(benchmark_file_dict)

    with PdfPages(output) as pdf:
        pdf.savefig(fig_sens)
        pdf.savefig(fig_a_eff)
        pdf.savefig(fig_ang_res)
        pdf.savefig(fig_e_res)


def parse_args():
    parser = ArgumentParser()
    parser.add_argument("--input_irfs", required=True, nargs="+")
    parser.add_argument("--input_benchmarks", required=True, nargs="+")
    parser.add_argument("-o", "--output", required=True)
    return parser.parse_args()


def main_from_snakemake(snakemake):
    main(snakemake.input.irfs, snakemake.input.benchmarks, snakemake.output[0])


if "snakemake" in globals():
    main_from_snakemake(snakemake)  # noqa: F821
elif __name__ == "__main__":
    args = parse_args()
    main(args.input_irfs, args.input_benchmarks, args.output)
