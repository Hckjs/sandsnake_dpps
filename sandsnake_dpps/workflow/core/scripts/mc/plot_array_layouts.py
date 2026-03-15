import matplotlib
import matplotlib.pyplot as plt
from argparse import ArgumentParser

from ctapipe.instrument import SubarrayDescription
from ctapipe.visualization import ArrayDisplay

# from ..scriptutils.logging import setup_logging

if matplotlib.get_backend() == "pgf":
    from matplotlib.backends.backend_pgf import PdfPages
else:
    from matplotlib.backends.backend_pdf import PdfPages

# log = logging.getLogger(__name__)


def main(input_path, output):
    figs = []
    for file in input_path:
        fig, ax = plt.subplots(1, 1, constrained_layout=True)
        ax = ax.ravel()
        subarray = SubarrayDescription.from_hdf(file)
        disp = ArrayDisplay(subarray, axes=ax)
        disp.add_labels()
        figs.append(fig)

    with PdfPages(output) as pdf:
        for fig in figs:
            pdf.savefig(fig)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-i", "--input-path", nargs="+", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("--log-file")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()
    # setup_logging(logfile=args.log_file, verbose=args.verbose)

    main(args.input_path, args.output)
