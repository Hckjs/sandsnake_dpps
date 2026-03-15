import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from argparse import ArgumentParser

from ctapipe.io import TableLoader
from ctapipe.instrument import SubarrayDescription

# from ..scriptutils.logging import setup_logging
import dl1_data_plots as dl1_plots

if matplotlib.get_backend() == "pgf":
    from matplotlib.backends.backend_pgf import PdfPages
else:
    from matplotlib.backends.backend_pdf import PdfPages

# log = logging.getLogger(__name__)


def calc_y_for_reco_vs_true(reco, true):
    out = np.zeros_like(true, dtype=true.dtype)
    return np.divide((reco - true), true, where=true != 0, out=out)


def main(input_path, output):
    tel_lst = "LST_LST_LSTcam"
    tel_mst = "MST_MST_NectarCam"
    subarray = SubarrayDescription.from_hdf(input_path)
    geom_lst = subarray.tel[1].camera.geometry
    geom_mst = subarray.tel[5].camera.geometry

    loader = TableLoader(input_path)
    tel_tables = loader.read_telescope_events_by_type(
        telescopes=[tel_lst, tel_mst],
        start=None,
        stop=None,
        dl1_images=True,
        dl1_parameters=True,
        dl1_muons=False,
        dl2=True,
        simulated=True,
        true_images=True,
        true_parameters=True,
        instrument=True,
        observation_info=False,
    )

    subarray_table = loader.read_subarray_events(
        start=None,
        stop=None,
        dl2=True,
        simulated=True,
        observation_info=False,
    )

    fig_hist_1, ax_hist_1 = plt.subplots(2, 2)
    fig_hist_2, ax_hist_2 = plt.subplots(2, 2)
    fig_hist_3, ax_hist_3 = plt.subplots(2, 1)

    fig_e_dep_1, ax_e_dep_1 = plt.subplots(2, 2)
    fig_e_dep_2, ax_e_dep_2 = plt.subplots(2, 2)
    fig_e_dep_3, ax_e_dep_3 = plt.subplots(2, 2)
    fig_e_dep_4, ax_e_dep_4 = plt.subplots(2, 2)

    fig_true_1, ax_true_1 = plt.subplots(2, 2)
    fig_true_2, ax_true_2 = plt.subplots(2, 2)

    fig_cd, ax_cd = plt.subplots(2, 2)

    # Histograms
    dl1_plots.hist_total_intensity(
        [
            tel_tables[tel_lst]["hillas_intensity"],
            tel_tables[tel_mst]["hillas_intensity"],
        ],
        ax_hist_1[0, 0],
        label=["LST", "MST"],
    )
    dl1_plots.hist_tels_with_trigger(
        subarray_table["tels_with_trigger"].value.sum(axis=1),
        ax_hist_1[0, 1],
    )
    dl1_plots.hist_trigger_per_tel(
        [
            subarray_table["tels_with_trigger"].value[:, :4].sum(axis=1),
            subarray_table["tels_with_trigger"].value[:, 4:].sum(axis=1),
        ],
        ax_hist_1[1, 0],
        label=["LST", "MST"],
    )
    dl1_plots.hist_selected_pixels(
        [
            tel_tables[tel_lst]["image_mask"].value.sum(axis=1),
            tel_tables[tel_mst]["image_mask"].value.sum(axis=1),
        ],
        ax_hist_1[1, 1],
        label=["LST", "MST"],
    )
    ics_per_tel_type = {}
    for tel_type in [tel_lst, tel_mst]:
        true_images = tel_tables[tel_type]["true_image"].value
        image_masks = tel_tables[tel_type]["image_mask"].value
        cleaned_true_images = np.where(image_masks == 0, 0, true_images)
        ics_per_tel_type[tel_type] = np.divide(
            cleaned_true_images.sum(axis=1),
            tel_tables[tel_type]["true_image_sum"].value,
        )

    dl1_plots.hist_identified_cherenkov_signal(
        ics_per_tel_type.values(),
        ax_hist_2[0, 0],
        label=["LST", "MST"],
    )
    dl1_plots.hist_tel_impact_distance(
        [
            tel_tables[tel_lst]["HillasReconstructor_tel_impact_distance"],
            tel_tables[tel_mst]["HillasReconstructor_tel_impact_distance"],
        ],
        ax_hist_2[0, 1],
        label=["LST", "MST"],
    )
    dl1_plots.hist_hillas_psi(
        [
            tel_tables[tel_lst]["hillas_psi"],
            tel_tables[tel_mst]["hillas_psi"],
        ],
        ax_hist_2[1, 0],
        label=["LST", "MST"],
    )
    dl1_plots.hist_leakage(
        [
            tel_tables[tel_lst]["leakage_intensity_width_2"],
            tel_tables[tel_mst]["leakage_intensity_width_2"],
        ],
        ax_hist_2[1, 1],
        label=["LST", "MST"],
    )
    dl1_plots.hist_concentration_core(
        [
            tel_tables[tel_lst]["concentration_core"],
            tel_tables[tel_mst]["concentration_core"],
        ],
        ax_hist_3[0],
        label=["LST", "MST"],
    )
    n_islands_lst = tel_tables[tel_lst]["morphology_n_islands"]
    n_islands_mst = tel_tables[tel_mst]["morphology_n_islands"]
    dl1_plots.hist_n_islands(
        [
            np.where(n_islands_lst < 0, 0, n_islands_lst),
            np.where(n_islands_mst < 0, 0, n_islands_mst),
        ],
        ax_hist_3[1],
        label=["LST", "MST"],
    )

    # Energy dependet plots
    tel_dict = {
        "LST": tel_lst,
        "MST": tel_mst,
    }

    for i, (tel_str, tel) in enumerate(tel_dict.items()):
        dl1_plots.hillas_intensity_per_energy(
            tel_tables[tel]["true_energy"],
            tel_tables[tel]["hillas_intensity"],
            fig_e_dep_1,
            ax_e_dep_1[0, i],
            tel_str,
        )

        dl1_plots.ics_per_energy(
            tel_tables[tel]["true_energy"],
            ics_per_tel_type[tel],
            fig_e_dep_1,
            ax_e_dep_1[1, i],
            tel_str,
        )

        dl1_plots.hillas_wl_per_energy(
            tel_tables[tel]["true_energy"],
            tel_tables[tel]["hillas_width"].value
            / tel_tables[tel]["hillas_length"].value,
            fig_e_dep_1,
            ax_e_dep_2[0, i],
            tel_str,
        )

        dl1_plots.hillas_wl_true_per_energy(
            tel_tables[tel]["true_energy"],
            tel_tables[tel]["true_hillas_width"].value
            / tel_tables[tel]["true_hillas_length"].value,
            fig_e_dep_1,
            ax_e_dep_2[1, i],
            tel_str,
        )

        dl1_plots.hillas_psi_per_energy(
            tel_tables[tel]["true_energy"],
            tel_tables[tel]["hillas_psi"],
            fig_e_dep_1,
            ax_e_dep_3[0, i],
            tel_str,
        )

        dl1_plots.intensity_max_per_energy(
            tel_tables[tel]["true_energy"],
            tel_tables[tel]["intensity_max"],
            fig_e_dep_1,
            ax_e_dep_3[1, i],
            tel_str,
        )

        dl1_plots.peak_time_max_min_per_energy(
            tel_tables[tel]["true_energy"],
            tel_tables[tel]["peak_time_max"].value
            - tel_tables[tel]["peak_time_min"].value,
            fig_e_dep_1,
            ax_e_dep_4[0, i],
            tel_str,
        )

        dl1_plots.tel_impact_distance_per_energy(
            tel_tables[tel]["true_energy"],
            tel_tables[tel]["HillasReconstructor_tel_impact_distance"],
            fig_e_dep_1,
            ax_e_dep_4[1, i],
            tel_str,
        )

        ## Reco vs. True ##
        dl1_plots.hillas_intensity_true(
            tel_tables[tel]["true_hillas_intensity"],
            calc_y_for_reco_vs_true(
                tel_tables[tel]["hillas_intensity"],
                tel_tables[tel]["true_hillas_intensity"],
            ),
            fig_true_1,
            ax_true_1[0, i],
            tel_str,
        )

        dl1_plots.concentration_core_true(
            tel_tables[tel]["true_concentration_core"],
            calc_y_for_reco_vs_true(
                tel_tables[tel]["concentration_core"],
                tel_tables[tel]["true_concentration_core"],
            ),
            fig_true_1,
            ax_true_1[1, i],
            tel_str,
        )

        n_islands = tel_tables[tel]["morphology_n_islands"]
        dl1_plots.n_islands_true(
            tel_tables[tel]["true_morphology_n_islands"],
            calc_y_for_reco_vs_true(
                np.float64(np.where(n_islands < 0, 0, n_islands)),
                np.float64(tel_tables[tel]["true_morphology_n_islands"]),
            ),
            fig_true_2,
            ax_true_2[0, i],
            tel_str,
        )

        dl1_plots.tel_impact_distance_true(
            tel_tables[tel]["true_impact_distance"],
            calc_y_for_reco_vs_true(
                tel_tables[tel]["HillasReconstructor_tel_impact_distance"],
                tel_tables[tel]["true_impact_distance"],
            ),
            fig_true_2,
            ax_true_2[1, i],
            tel_str,
        )

    ## Camera Display ##
    # Might be better to plot per tel_id and not per telescope type
    dl1_plots.cam_intensity_mean(
        geom_lst,
        tel_tables[tel_lst]["image"].value.mean(axis=0),
        ax_cd[0, 0],
    )

    dl1_plots.cam_intensity_mean(
        geom_mst,
        tel_tables[tel_mst]["image"].value.mean(axis=0),
        ax_cd[0, 1],
    )

    dl1_plots.cam_intensity_std(
        geom_lst,
        tel_tables[tel_lst]["image"].value.std(axis=0),
        ax_cd[1, 0],
    )

    dl1_plots.cam_intensity_std(
        geom_mst,
        tel_tables[tel_mst]["image"].value.std(axis=0),
        ax_cd[1, 1],
    )

    fig_hist_1.tight_layout()
    fig_hist_2.tight_layout()
    fig_hist_3.tight_layout()
    fig_cd.tight_layout()
    fig_e_dep_1.tight_layout()
    fig_e_dep_2.tight_layout()
    fig_e_dep_3.tight_layout()
    fig_e_dep_4.tight_layout()
    fig_true_1.tight_layout()
    fig_true_2.tight_layout()

    with PdfPages(output) as pdf:
        pdf.savefig(fig_hist_1)
        pdf.savefig(fig_hist_2)
        pdf.savefig(fig_hist_3)
        pdf.savefig(fig_e_dep_1)
        pdf.savefig(fig_e_dep_2)
        pdf.savefig(fig_e_dep_3)
        pdf.savefig(fig_e_dep_4)
        pdf.savefig(fig_true_1)
        pdf.savefig(fig_true_2)
        pdf.savefig(fig_cd)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-i", "--input-path", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("--log-file")
    args = parser.parse_args()
    # setup_logging(logfile=args.log_file, verbose=args.verbose)

    main(args.input_path, args.output)
