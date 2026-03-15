from pathlib import Path

# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------


def _ext(input_key: str, re_none: bool = True):
    val = config.get("inputs", {}).get(input_key)
    if val in (None, "", "null"):
        if re_none:
            return None
        return ""

    return str(Path(val).expanduser().resolve())


def bind_wildcards(provider, **overrides):
    def _call(wc):
        data = dict(wc)
        data.update(overrides)
        if isinstance(provider, str):
            return provider.format(**data)
        path = provider(wc)
        return path.format(**data)

    return _call


def _select_input(input_key: str) -> str:
    """
    Return "stage" if external input is given for config_key,
    else return "core".
    """

    ext_path = _ext(input_key)
    if ext_path:
        return "stage"

    return "core"


# ------------------------------------------------------------------
# MC providers
# ------------------------------------------------------------------


def mc_simtel_provider(wc):
    # DL1 reprocessing
    if _ext("mc_dl1") and "mc_dl1" in config.get("targets", []):
        return mc_dl1_provider(wc)

    base = _select_input("mc_simtel")
    return TARGETS_MC_SIMTEL(base, resolve=False)


def mc_dl1_provider(wc):
    base = _select_input("mc_dl1")
    return TARGETS_MC_DL1(base, resolve=False)


def mc_dl1_merged_provider(wc):
    base = _select_input("mc_dl1_merged")
    return TARGETS_MC_DL1_MERGED(base, resolve=False)


def mc_dl2_provider(wc):
    base = _select_input("mc_dl2")
    return TARGETS_MC_DL2(base, resolve=False)


# ------------------------------------------------------------------
# RF model provider
# ------------------------------------------------------------------


def rf_energy_regressor_provider(wc):
    base = _select_input("rf_energy_regressor")
    return TARGETS_RF_ENERGY_REGRESSOR(base, resolve=False)


def rf_particle_classifier_provider(wc):
    base = _select_input("rf_particle_classifier")
    return TARGETS_RF_PARTICLE_CLASSIFIER(base, resolve=False)


def rf_geometry_reconstructor_provider(wc):
    base = _select_input("rf_geometry_reconstructor")
    return TARGETS_RF_GEOMETRY_RECONSTRUCTOR(base, resolve=False)


# ------------------------------------------------------------------
# IRFs provider
# ------------------------------------------------------------------


def irfs_provider(wc):
    base = _select_input("irfs")
    return TARGETS_IRFS(base, resolve=False)


def benchmarks_provider(wc):
    base = _select_input("irfs")
    return TARGETS_BENCHMARKS(base, resolve=False)


def cuts_provider(wc):
    base = _select_input("irfs")
    return TARGETS_CUTS(base, resolve=False)


# ------------------------------------------------------------------
# Internal provider
# ------------------------------------------------------------------


def mc_dl1_split_provider(wc):
    train_size = config.get("train_size")
    cuts_size = config.get("cuts_size")
    split = getattr(wc, "split", None)
    particle = getattr(wc, "particle", None)

    base_dir = _ext("mc_dl1")
    # if Dl1 reprocessing
    if base_dir and "mc_dl1" in config.get("targets", []):
        base_dir = PATHS["core:mc_dl1"]

    parent_dir = f"/zen_{wc.zen}/az_{wc.az}/{particle}"
    filenames = get_filenames(wc.zen, wc.az, particle)

    split_path = Path(base_dir + parent_dir)
    split_files = [split_path / f"{fname}.dl1.h5" for fname in filenames]

    num_files = len(split_files)
    train_files = split_files[: int(num_files * train_size)]
    test_files = split_files[int(num_files * train_size) :]

    # TODO: Make a check that splits are not empty

    if particle == "gamma_diffuse":
        if split == "train_en":
            return train_files[: int(len(train_files) * 0.4)]
        if split == "train_cl_disp":
            return train_files[int(len(train_files) * 0.4) :]
        if split == "test_cuts":
            return test_files[: int(len(test_files) * cuts_size)]
        if split == "test_irfs":
            return test_files[int(len(test_files) * cuts_size) :]

    if particle == "proton":
        if split == "train_cl_disp":
            return train_files
        if split == "test_cuts":
            return test_files[: int(len(test_files) * cuts_size)]
        if split == "test_irfs":
            return test_files[int(len(test_files) * cuts_size) :]

    if particle == "gamma":
        if split == "test_cuts":
            return split_files[: int(num_files * cuts_size)]
        if split == "test_irfs":
            return split_files[int(num_files * cuts_size) :]

    if particle == "electron":
        if split == "test_cuts":
            return split_files[: int(num_files * cuts_size)]
        if split == "test_irfs":
            return split_files[int(num_files * cuts_size) :]
