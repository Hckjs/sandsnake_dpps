import matplotlib
import matplotlib.pyplot as plt
from argparse import ArgumentParser
import numpy as np

from ctapipe.io import TableLoader
import dl2_rf_performance_plots as dl2_plots

# from ..scriptutils.logging import setup_logging

if matplotlib.get_backend() == "pgf":
    from matplotlib.backends.backend_pgf import PdfPages
else:
    from matplotlib.backends.backend_pdf import PdfPages

# log = logging.getLogger(__name__)


def main(
    input_gammas, input_protons, disp_model, e_reg_model, clf_model, config_path, output
):
    fig_dir_reco, ax_dir_reco = plt.subplots(2, 2)
    fig_e_reco, ax_e_reco = plt.subplots(2, 2)
    fig_gh, ax_gh = plt.subplots(2, 2)
    # fig_feat, ax_feat = plt.subplots(1, 2)

    with TableLoader(input_gammas) as loader_gammas:
        gammas = loader_gammas.read_subarray_events(
            start=None,
            stop=None,
            dl2=True,
            simulated=True,
            observation_info=True,
        )
        nans_gammas = np.isnan(gammas["RandomForestRegressor_energy"].value)
        gammas = gammas[~nans_gammas]

        dl2_plots.direction_results(gammas, ax_dir_reco)
        dl2_plots.plot_energy_results(gammas, ax_e_reco)
        gamma_cl_prediction = gammas["RandomForestClassifier_prediction"].value
        gamma_true_energy = gammas["true_energy"].value
        del gammas
        del nans_gammas

    fig_dir_reco.tight_layout()
    fig_e_reco.tight_layout()

    with TableLoader(input_protons) as loader_protons:
        protons = loader_protons.read_subarray_events(
            start=None,
            stop=None,
            dl2=True,
            simulated=True,
            observation_info=False,
        )

        nans_protons = np.isnan(protons["RandomForestRegressor_energy"].value)
        protons = protons[~nans_protons]
        protons_cl_prediction = protons["RandomForestClassifier_prediction"].value
        protons_true_energy = protons["true_energy"].value
        del protons
        del nans_protons

    dl2_plots.plot_roc_gamma(
        gamma_cl_prediction,
        gamma_true_energy,
        protons_cl_prediction,
        protons_true_energy,
        ax_gh,
    )

    fig_gh.tight_layout()

    figs_imp = dl2_plots.plot_models_feature_importances(
        e_reg_model,
        disp_model,
        clf_model,
        config_path,
    )
    # dl2_plots.plot_features(gammas, protons, ax_feat)
    # fig_feat.tight_layout()

    with PdfPages(output) as pdf:
        pdf.savefig(fig_gh)
        pdf.savefig(fig_dir_reco)
        pdf.savefig(fig_e_reco)
        pdf.savefig(figs_imp)
        # pdf.savefig(fig_feat)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-g", "--gammas", required=True)
    parser.add_argument("-p", "--protons", required=True)
    parser.add_argument("-d", "--disp-model", required=True)
    parser.add_argument("-e", "--e-reg-model", required=True)
    parser.add_argument("-clf", "--clf-model", required=True)
    parser.add_argument("-c", "--config", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("--log-file")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()
    # setup_logging(logfile=args.log_file, verbose=args.verbose)

    main(
        args.gammas,
        args.protons,
        args.disp_model,
        args.e_reg_model,
        args.clf_model,
        args.config,
        args.output,
    )
