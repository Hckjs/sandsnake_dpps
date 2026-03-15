# Stereo-combiner plugin paths
from pathlib import Path

STEREO_DIR = BUILD_DIR / "stereo_combiner" / "{run_tag}"

paths_update(
    "stereo",
    {
        # intermediate / primary products
        "combined_dl2": STEREO_DIR / "combined_dl2.h5",
        # products often used as "final"
        "irfs_combined": STEREO_DIR / "irfs_combined.fits.gz",
        "cuts_table": STEREO_DIR / "optimized_cuts.ecsv",
        # optional plots
        "sens_plot": STEREO_DIR / "sensitivity.png",
    },
)

stereo_combiners = [
    "StereoMeanCombiner",
    "StereoDispCombiner",
    "StereoDispCombinerSS",
    "StereoKMeansCombiner",
    "StereoDBScanCombiner",
]


def DL2_STEREO_COMB_DATA(wildcards):
    combiner_path = OUTDIRS["stereo_combiner"]
    file_list = [
        combiner_path / f"mc/zen_20/az_0/gamma/test_irfs/"
        f"gamma_zen_20_az_0_test_irfs_ac_full_array_{combiner}.dl2.h5"
        for combiner in stereo_combiners
    ]
    return file_list


def IRFS_STEREO_COMB_DATA(file_type):
    combiner_path = OUTDIRS["stereo_combiner"]
    file_list = [
        combiner_path
        / f"mc/zen_20/az_0/irfs/{file_type}_zen_20_az_0_ac_full_array_{combiner}.fits.gz"
        for combiner in stereo_combiners
    ]
    return file_list
