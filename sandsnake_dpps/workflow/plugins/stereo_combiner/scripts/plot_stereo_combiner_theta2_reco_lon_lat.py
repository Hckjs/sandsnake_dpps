from argparse import ArgumentParser
import matplotlib
from ctapipe.io import TableLoader
from pathlib import Path

from .stereo_combiner_plots import (
    plot_theta2,
    plot_lon_lat,
    stack_theta_hist,
    stack_lon_lat_hist,
)

if matplotlib.get_backend() == "pgf":
    from matplotlib.backends.backend_pgf import PdfPages
else:
    from matplotlib.backends.backend_pdf import PdfPages


def main(dl2_gammas, output):
    theta_hist_dict = {}
    lon_lat_hist_dict = {}

    # loop over combiners
    for gamma_file in dl2_gammas:
        combiner = Path(gamma_file).stem.split("_")[-1].split(".")[0]
        with TableLoader(input_url=gamma_file) as loader:
            chunk_iterator = loader.read_subarray_events_chunked(
                chunk_size=1000000,
                dl2=True,
                simulated=True,
                observation_info=True,
            )

            hist_theta = None
            hist_lon_lat = None
            for _, _, dl2_table in chunk_iterator:
                # Mask of tel_ids
                multiplicity_mask = (dl2_table["disp_telescopes"].value).sum(-1) > 1
                dl2_table = dl2_table[multiplicity_mask]

                theta_hist_chunk, theta_hist_edges = stack_theta_hist(dl2_table)
                if hist_theta is None:
                    hist_theta = theta_hist_chunk
                else:
                    hist_theta += theta_hist_chunk

                lon_lat_hist_chunk, extent = stack_lon_lat_hist(dl2_table)
                if hist_lon_lat is None:
                    hist_lon_lat = lon_lat_hist_chunk
                else:
                    hist_lon_lat += lon_lat_hist_chunk

            theta_hist_dict[combiner] = (hist_theta, theta_hist_edges)
            lon_lat_hist_dict[combiner] = hist_lon_lat

    # Create Figs
    with PdfPages(output) as pdf:
        fig_theta = plot_theta2(theta_hist_dict)
        fig_lon_lat = plot_lon_lat(lon_lat_hist_dict, extent)

        pdf.savefig(fig_theta)
        pdf.savefig(fig_lon_lat)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--input", required=True, nargs="+")
    parser.add_argument("-o", "--output", required=True)
    args = parser.parse_args()

    main(args.input, args.output)
