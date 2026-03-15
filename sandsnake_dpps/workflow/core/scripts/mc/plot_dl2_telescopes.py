from argparse import ArgumentParser
import matplotlib
from ctapipe.io import TableLoader

from dl2_plots import (
    plot_theta_2_tel,
    plot_reco_alt_az_tel,
    stack_theta_hist,
    stack_alt_az_hist,
)

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
        chunk_iterator = loader.read_telescope_events_by_type_chunked(
            chunk_size=1000000,
            dl2=True,
            simulated=True,
            observation_info=True,
            dl1_parameters=False,
            dl1_images=False,
            true_images=False,
        )

        hist_theta_total = [None, None]
        hist_alt_az_total = [None, None]
        tel_types = [None, None]
        for chunk, (_, _, tel_tables) in enumerate(chunk_iterator):
            for i, (tel_type, dl2_table) in enumerate(tel_tables.items()):
                # mask = dl2_table["true_energy"].value > 2
                # print(len(dl2_table), mask.sum())
                # dl2_table = dl2_table[mask]
                breakpoint()
                theta_hist = stack_theta_hist(dl2_table)
                if hist_theta_total[i] is None:
                    hist_theta_total[i] = theta_hist
                else:
                    hist_theta_total[i] += theta_hist

                alt_az_hist = stack_alt_az_hist(dl2_table)
                if hist_alt_az_total[i] is None:
                    hist_alt_az_total[i] = alt_az_hist
                else:
                    hist_alt_az_total[i] += alt_az_hist

                tel_types[i] = tel_type

            print(chunk)
            if chunk == 0:
                break

        with PdfPages(output) as pdf:
            for i, tel_type in enumerate(tel_types):
                fig_theta_2 = plot_theta_2_tel(hist_theta_total[i], tel_type)
                fig_reco_alt_az = plot_reco_alt_az_tel(hist_alt_az_total[i], tel_type)

                pdf.savefig(fig_theta_2)
                pdf.savefig(fig_reco_alt_az)


if __name__ == "__main__":
    main(**vars(args))
