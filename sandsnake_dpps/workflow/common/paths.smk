# ----------------------------------------------------------------------
# Repo / build roots
# ----------------------------------------------------------------------

WORKFLOW_DIR = Path(workflow.basedir).expanduser().resolve()
BUILD_DIR = Path(config.get("build_dir", "build")).expanduser().resolve()

CORE_DIR = WORKFLOW_DIR / "core"
PLUGINS_DIR = WORKFLOW_DIR / "plugins"

CORE_SCRIPTS_DIR = CORE_DIR / "scripts"
CORE_ENVS_DIR = CORE_DIR / "envs"
CORE_CONFIGS_DIR = CORE_DIR / "configs"

USER_CONFIGS_DIR = config.get("user_config_dir", None)
USER_ENVS_DIR = config.get("user_env_dir", None)

RUN_TAG = config.get("run_tag", "default")

os.environ["PYTHONPATH"] = str(WORKFLOW_DIR)

# ----------------------------------------------------------------------
# Directory roots
# ----------------------------------------------------------------------

RUN_DIR = str(BUILD_DIR / RUN_TAG)
OUTDIRS = {
    "dl1": f"{RUN_DIR}/dl1",
    "dl2": f"{RUN_DIR}/dl2",
    "dl3": f"{RUN_DIR}/dl3",
    "mc": f"{RUN_DIR}/mc",
    "models": f"{RUN_DIR}/models",
    "irfs": f"{RUN_DIR}/irfs",
    "plots": f"{RUN_DIR}/plots",
    "staging": f"{RUN_DIR}/staging",
    "plugins": f"{RUN_DIR}/plugins",
    "configs": f"{RUN_DIR}/configs",
    "envs": f"{RUN_DIR}/envs",
}

# ----------------------------------------------------------------------
# Paths
# ----------------------------------------------------------------------

PATHS = {
    # ------------------
    # Obs Dirs
    # ------------------
    "core:obs_dl1": OUTDIRS["dl1"],
    "core:obs_dl2": OUTDIRS["dl2"],
    "core:obs_dl3": OUTDIRS["dl3"],
    "stage:obs_dl1": f'{OUTDIRS["staging"]}/dl1',
    "stage:obs_dl2": f'{OUTDIRS["staging"]}/dl2',
    "stage:obs_dl3": f'{OUTDIRS["staging"]}/dl3',
    # ------------------
    # MC Dirs
    # ------------------
    "core:mc_simtel": f'{OUTDIRS["mc"]}/simtel',
    "core:mc_dl1": f'{OUTDIRS["mc"]}/dl1',
    "core:mc_dl1_merged": f'{OUTDIRS["mc"]}/dl1_merged',
    "core:mc_dl2": f'{OUTDIRS["mc"]}/dl2',
    "stage:mc_simtel": f'{OUTDIRS["staging"]}/mc/simtel',
    "stage:mc_dl1": f'{OUTDIRS["staging"]}/mc/dl1',
    "stage:mc_dl1_splits": f'{OUTDIRS["staging"]}/mc/dl1_splits',
    "stage:mc_dl1_merged": f'{OUTDIRS["staging"]}/mc/dl1_merged',
    "stage:mc_dl2": f'{OUTDIRS["staging"]}/mc/dl2',
    # ------------------
    # IRFs
    # ------------------
    "core:irfs": OUTDIRS["irfs"],
    "stage:irfs": f'{OUTDIRS["staging"]}/irfs',
    # ------------------
    # RF Models
    # ------------------
    "core:models": OUTDIRS["models"],
    "stage:models": f'{OUTDIRS["staging"]}/models',
    # ------------------
    # Plots
    # ------------------
    "core:plots": OUTDIRS["plots"],
    # ------------------
    # Configs
    # ------------------
    "core:config:merge": f'{OUTDIRS["configs"]}/core/merge.yml',
    "core:config:process": f'{OUTDIRS["configs"]}/core/process.yml',
    "core:config:train_models": f'{OUTDIRS["configs"]}/core/train_models.yml',
    "core:config:apply_models": f'{OUTDIRS["configs"]}/core/apply_models.yml',
    "core:config:optimize_cuts": f'{OUTDIRS["configs"]}/core/optimize_cuts.yml',
    "core:config:compute_irf": f'{OUTDIRS["configs"]}/core/compute_irf.yml',
    # ------------------
    # Conda Environments
    # ------------------
    "core:env:ctapipe": f'{OUTDIRS["envs"]}/core/ctapipe.yaml',
    "core:env:plotting": f'{OUTDIRS["envs"]}/core/plotting.yaml',
    "core:env:gammapy": f'{OUTDIRS["envs"]}/core/gammapy.yaml',
    # ------------------
    # Templates (suffixes only; no RUN_DIR prefix)
    # ------------------
    # Main
    "core:template:mc_dl1": "/zen_{zen}/az_{az}/{particle}",
    "core:template:particle_split": "/zen_{zen}/az_{az}/{particle}/{split}",
    "core:template:mc_dl1_merged": "/zen_{zen}/az_{az}/{particle}/{split}/{particle}_zen_{zen}_az_{az}_{split}_merged.dl1.h5",
    "core:template:rf_energy_regressor": "/zen_{zen}/az_{az}/rf_energy_regressor_zen_{zen}_az_{az}.pkl",
    "core:template:rf_particle_classifier": "/zen_{zen}/az_{az}/rf_particle_classifier_zen_{zen}_az_{az}.pkl",
    "core:template:rf_geometry_reconstructor": "/zen_{zen}/az_{az}/rf_geometry_reconstructor_zen_{zen}_az_{az}.pkl",
    "core:template:mc_dl2": "/zen_{zen}/az_{az}/{particle}/{split}/{particle}_zen_{zen}_az_{az}_{split}.dl2.h5",
    "core:template:rf_performance_plots": "/zen_{zen}/az_{az}/model_performance_zen_{zen}_az_{az}.pdf",
    "core:template:cuts": "/zen_{zen}/az_{az}/cuts_zen_{zen}_az_{az}_obs_{obstime}_hours.fits",
    "core:template:irfs": "/zen_{zen}/az_{az}/irfs_zen_{zen}_az_{az}_obs_{obstime}_hours.fits.gz",
    "core:template:benchmarks": "/zen_{zen}/az_{az}/benchmarks_zen_{zen}_az_{az}_obs_{obstime}_hours.fits.gz",
    "core:template:irf_plots": "/zen_{zen}/az_{az}/irfs_zen_{zen}_az_{az}_obs_{obstime}_hours.pdf",
    # Misc
    "core:template:mc_dl1_applied_e_reg": "/zen_{zen}/az_{az}/{particle}/train_cl_disp/"
    "{particle}_zen_{zen}_az_{az}_train_cl_disp_applied.dl1.h5",
}

OUTPATHS = {
    "mc_simtel": PATHS["core:mc_simtel"]
    + PATHS["core:template:mc_dl1"]
    + "/{filename}.simtel.gz",
    "mc_dl1": PATHS["core:mc_dl1"]
    + PATHS["core:template:mc_dl1"]
    + "/{filename}.dl1.h5",
    "mc_dl1_split": PATHS["stage:mc_dl1"] + PATHS["core:template:particle_split"],
    "mc_dl1_merged": PATHS["core:mc_dl1_merged"] + PATHS["core:template:mc_dl1_merged"],
    "apply_energy_regressor_train": PATHS["core:mc_dl1"]
    + PATHS["core:template:mc_dl1_applied_e_reg"],
    "train_energy_regressor": PATHS["core:models"]
    + PATHS["core:template:rf_energy_regressor"],
    "train_disp_reconstructor": PATHS["core:models"]
    + PATHS["core:template:rf_geometry_reconstructor"],
    "train_particle_classifier": PATHS["core:models"]
    + PATHS["core:template:rf_particle_classifier"],
    "mc_dl2": PATHS["core:mc_dl2"] + PATHS["core:template:mc_dl2"],
    "rf_performance_plots": PATHS["core:plots"]
    + PATHS["core:template:rf_performance_plots"],
    "cuts": PATHS["core:irfs"] + PATHS["core:template:cuts"],
    "irfs": PATHS["core:irfs"] + PATHS["core:template:irfs"],
    "benchmarks": PATHS["core:irfs"] + PATHS["core:template:benchmarks"],
    "plot_irfs": PATHS["core:plots"] + PATHS["core:template:irf_plots"],
}
# ----------------------------------------------------------------------
# Path resolution helpers
# ----------------------------------------------------------------------


def paths_update(namespace: str, mapping: dict):
    """
    Extend PATHS under a plugin namespace.

    Example:
        paths_update("fermi", {"sensitivity": Path(...)})
    """
    if not namespace or ":" in namespace:
        raise ValueError("Invalid namespace for PATHS extension")

    for k, v in mapping.items():
        kk = f"{namespace}:{k}"
        if kk in PATHS:
            raise KeyError(f"PATHS key collision: {kk}")
        if not isinstance(v, Path):
            v = Path(v)
        PATHS[kk] = v


def log_path(output_path: str, suffix: str):
    p = Path(output_path)

    # template file
    if p.suffix:
        base = p.parent
        filename = p.stem
    # template dir
    else:
        base = p
        filename = p.name
    suf = suffix if suffix.startswith(".") else f".{suffix}"

    return base / "logs" / f"{filename}{suf}"
