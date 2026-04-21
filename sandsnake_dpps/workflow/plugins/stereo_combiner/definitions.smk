STEREO_COMBINER_DIR = PLUGINS_DIR / "stereo_combiner"
STEREO_COMBINER_SCRIPTS_DIR = STEREO_COMBINER_DIR / "scripts"

STEREO_COMBINER_OUTDIR = f'{OUTDIRS["plugins"]}/stero_combiner'
STEREO_COMBINER_PATHS = {
    "mc_dl2": STEREO_COMBINER_OUTDIR
    + "/mc/dl2/zen_20/az_180/{particle}"
    + "/{split}/{particle}_zen_20_az_180_{split}_comb_{combiner}.dl2.h5",
    "cuts": STEREO_COMBINER_OUTDIR
    + "/irfs/zen_20/az_180"
    + "/cuts_zen_20_az_180_obs_50_hours_comb_{combiner}.fits",
    "irfs": STEREO_COMBINER_OUTDIR
    + "/irfs/zen_20/az_180"
    + "/irfs_zen_20_az_180_obs_50_hours_comb_{combiner}.fits",
    "benchmarks": STEREO_COMBINER_OUTDIR
    + "/irfs/zen_20/az_180"
    + "/benchmarks_zen_20_az_180_obs_50_hours_comb_{combiner}.fits",
    "plot_reco_lon_lat": STEREO_COMBINER_OUTDIR
    + "/plots/stereo_combiner_theta2_reco_lon_lat.pdf",
    "plot_irfs": STEREO_COMBINER_OUTDIR + "/plots/stereo_combiner_irfs.pdf",
}

paths_update("stereo_combiner", STEREO_COMBINER_PATHS)

stereo_combiners = [
    "StereoMeanCombiner",
    "StereoDispCombiner",
    "StereoKMeansCombiner",
    "StereoDBScanCombiner",
]


def STEREO_COMBINER_MC_DL2(wildcards):
    raw_path = PATHS["stereo_combiner:mc_dl2"]
    file_list = [
        raw_path.format(particle="gamma", split="test_irfs", combiner=combiner)
        for combiner in stereo_combiners
    ]
    return file_list


def STEREO_COMBINER_IRFS(file_type):
    raw_path = PATHS[f"stereo_combiner:{file_type}"]
    file_list = [raw_path.format(combiner=combiner) for combiner in stereo_combiners]
    return file_list
