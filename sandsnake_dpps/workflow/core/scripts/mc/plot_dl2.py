from argparse import ArgumentParser
import matplotlib
from ctapipe.io import TableLoader

from dl2_plots import plot_theta_2, plot_reco_lon_lat

if matplotlib.get_backend() == "pgf":
    from matplotlib.backends.backend_pgf import PdfPages
else:
    from matplotlib.backends.backend_pdf import PdfPages


parser = ArgumentParser()
parser.add_argument("--dl2-gammas", required=True)
parser.add_argument("-o", "--output", required=True)
args = parser.parse_args()


def main(dl2_gammas, output):
    with TableLoader(input_url=dl2_gammas) as loader:
        dl2_table = loader.read_subarray_events(
            dl2=True, simulated=True, observation_info=True
        )
        fig_theta_2 = plot_theta_2(dl2_table)
        fig_reco_alt_az = plot_reco_lon_lat(dl2_table)

        with PdfPages(output) as pdf:
            pdf.savefig(fig_theta_2)
            pdf.savefig(fig_reco_alt_az)


if __name__ == "__main__":
    main(**vars(args))
