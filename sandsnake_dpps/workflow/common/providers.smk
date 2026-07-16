# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------


def _ext(input_key: str, re_none: bool = True):
    inputs = config.get("inputs", None)
    if not isinstance(inputs, dict):
        raise ValueError("config['inputs'] must be a mapping of input names to paths")

    val = inputs.get(input_key)
    if val in (None, "", "null"):
        if re_none:
            return None
        return ""

    return str(Path(val).expanduser().resolve())


def _target(target_key: str) -> bool:
    targets = config.get("targets", None)
    if not isinstance(targets, dict):
        raise ValueError(
            "config['targets'] must be a mapping of target names to booleans"
        )
    return targets.get(target_key)


def _select_input(input_key: str) -> str:
    """
    Return "stage" if external input is given for config_key,
    else return "core".
    """

    ext_path = _ext(input_key)
    if ext_path:
        return "stage"

    return "core"


def bind_wildcards(provider, **overrides):
    def _call(wc):
        data = dict(wc)
        data.update(overrides)
        if isinstance(provider, str):
            return provider.format(**data)
        path = provider(wc)
        return path.format(**data)

    return _call


# ------------------------------------------------------------------
# MC providers
# ------------------------------------------------------------------


def mc_simtel_provider(wc):
    # DL1 reprocessing
    if _ext("mc_dl1") and _target("mc_dl1"):
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
    base = _select_input("models")
    return TARGETS_RF_ENERGY_REGRESSOR(base, resolve=False)


def rf_particle_classifier_provider(wc):
    base = _select_input("models")
    return TARGETS_RF_PARTICLE_CLASSIFIER(base, resolve=False)


def rf_geometry_reconstructor_provider(wc):
    base = _select_input("models")
    return TARGETS_RF_GEOMETRY_RECONSTRUCTOR(base, resolve=False)


# ------------------------------------------------------------------
# IRFs provider
# ------------------------------------------------------------------


def irfs_gamma_particle() -> str:
    """Return the configured gamma MC sample for cut optimization and IRFs."""
    irfs_config = config.get("irfs", {})
    if not isinstance(irfs_config, dict):
        raise ValueError("config['irfs'] must be a mapping when configured")

    gamma_particle = irfs_config.get("gamma_particle", "gamma")
    valid_particles = {"gamma", "gamma_diffuse"}
    if not isinstance(gamma_particle, str) or gamma_particle not in valid_particles:
        valid_values = ", ".join(sorted(valid_particles))
        raise ValueError(
            "config['irfs']['gamma_particle'] must be one of "
            f"{valid_values}; got {gamma_particle!r}"
        )

    return gamma_particle


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
    if base_dir and _target("mc_dl1"):
        base_dir = PATHS["core:mc_dl1"]

    parent_dir = f"/zen_{wc.zen}/az_{wc.az}/{particle}"
    filenames = get_filenames(wc.zen, wc.az, particle)

    split_path = Path(base_dir + parent_dir)
    split_files = [split_path / f"{fname}.dl1.h5" for fname in filenames]

    num_files = len(split_files)
    train_files = split_files[: int(num_files * train_size)]
    test_files = split_files[int(num_files * train_size) :]

    def non_empty(files):
        if not files:
            raise ValueError(
                f"MC DL1 split '{split}' for particle '{particle}' is empty."
            )
        return files

    if particle == "gamma_diffuse":
        if split == "train_en":
            return non_empty(train_files[: int(len(train_files) * 0.4)])
        if split == "train_cl_disp":
            return non_empty(train_files[int(len(train_files) * 0.4) :])
        if split == "test_cuts":
            return non_empty(test_files[: int(len(test_files) * cuts_size)])
        if split == "test_irfs":
            return non_empty(test_files[int(len(test_files) * cuts_size) :])

    if particle == "proton":
        if split == "train_cl_disp":
            return non_empty(train_files)
        if split == "test_cuts":
            return non_empty(test_files[: int(len(test_files) * cuts_size)])
        if split == "test_irfs":
            return non_empty(test_files[int(len(test_files) * cuts_size) :])

    if particle == "gamma":
        if split == "test_cuts":
            return non_empty(split_files[: int(num_files * cuts_size)])
        if split == "test_irfs":
            return non_empty(split_files[int(num_files * cuts_size) :])

    if particle == "electron":
        if split == "test_cuts":
            return non_empty(split_files[: int(num_files * cuts_size)])
        if split == "test_irfs":
            return non_empty(split_files[int(num_files * cuts_size) :])
