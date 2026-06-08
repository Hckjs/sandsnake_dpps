# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------


def _stage(input_path: str, output_path: str, mode: str = "symlink"):
    """
    Stage input_path to output_path.
    mode: "symlink" (default) or "copy"
    """
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)

    if out.exists() or out.is_symlink():
        out.unlink()

    if mode == "symlink":
        out.symlink_to(Path(input_path).resolve())
    elif mode == "copy":
        shutil.copy2(input_path, out)
    else:
        raise ValueError(f"Unknown staging mode: {mode}")


# ------------------------------------------------------------------
# MC staging
# ------------------------------------------------------------------
#
rule stage_mc_dl1:
    input:
        _ext("mc_dl1", re_none=False)
        + PATHS["core:template:mc_dl1"]
        + "/{filename}.DL1.h5",
    output:
        PATHS["stage:mc_dl1"] + PATHS["core:template:mc_dl1"] + "/{filename}.dl1.h5",
    run:
        _stage(input[0], output[0])


rule stage_mc_dl1_merged:
    input:
        _ext("mc_dl1_merged", re_none=False) + PATHS["core:template:mc_dl1_merged"],
    output:
        PATHS["stage:mc_dl1_merged"] + PATHS["core:template:mc_dl1_merged"],
    run:
        _stage(input[0], output[0])


rule stage_mc_dl2:
    input:
        _ext("mc_dl2", re_none=False) + PATHS["core:template:mc_dl2"],
    output:
        PATHS["stage:mc_dl2"] + PATHS["core:template:mc_dl2"],
    run:
        _stage(input[0], output[0])


# ------------------------------------------------------------------
# RF model staging
# ------------------------------------------------------------------


rule stage_rf_energy_regressor:
    input:
        _ext("models", re_none=False) + PATHS["core:template:rf_energy_regressor"],
    output:
        PATHS["stage:models"] + PATHS["core:template:rf_energy_regressor"],
    run:
        _stage(input[0], output[0])


rule stage_rf_particle_classifier:
    input:
        _ext("models", re_none=False) + PATHS["core:template:rf_particle_classifier"],
    output:
        PATHS["stage:models"] + PATHS["core:template:rf_particle_classifier"],
    run:
        _stage(input[0], output[0])


rule stage_rf_geometry_reconstructor:
    input:
        _ext("models", re_none=False) + PATHS["core:template:rf_geometry_reconstructor"],
    output:
        PATHS["stage:models"] + PATHS["core:template:rf_geometry_reconstructor"],
    run:
        _stage(input[0], output[0])


# ------------------------------------------------------------------
# IRF staging
# ------------------------------------------------------------------


rule stage_cuts:
    input:
        _ext("cuts", re_none=False) + PATHS["core:template:cuts"],
    output:
        PATHS["stage:irfs"] + PATHS["core:template:cuts"],
    run:
        _stage(input[0], output[0])


rule stage_benchmarks:
    input:
        _ext("benchmarks", re_none=False) + PATHS["core:template:benchmarks"],
    output:
        PATHS["stage:irfs"] + PATHS["core:template:benchmarks"],
    run:
        _stage(input[0], output[0])


rule stage_irfs:
    input:
        _ext("irfs", re_none=False) + PATHS["core:template:irfs"],
    output:
        PATHS["stage:irfs"] + PATHS["core:template:irfs"],
    run:
        _stage(input[0], output[0])


# ------------------------------------------------------------------
# Config/Env staging
# ------------------------------------------------------------------


def select_config(config_name: str, plugin: str) -> str:
    if USER_CONFIGS_DIR is not None:
        user_config_dir = Path(USER_CONFIGS_DIR).expanduser().resolve()
        user_config = user_config_dir / plugin / f"{config_name}.yml"

        if user_config.exists():
            return str(user_config)

        message = (
            f"Configured user_config_dir is set to {user_config_dir}, but no user "
            f"configuration file was found for plugin {plugin!r} and config "
            f"{config_name!r} at {user_config}. User configuration files must "
            "be placed at <user_config_dir>/<plugin>/<config_name>.yml, must "
            "use the .yml extension, and must be named exactly "
            f"{config_name}.yml. Expected path for this configuration: "
            f"{user_config}. Falling back to the default configuration file."
        )
        smk_logger.warning(message)

    if plugin == "core":
        DEFAULT_CONFIG_DIR = CORE_CONFIGS_DIR
    else:
        DEFAULT_CONFIG_DIR = PLUGINS_DIR / plugin / "configs"

    default_config = DEFAULT_CONFIG_DIR / f"{config_name}.yml"
    if not default_config.exists():
        raise FileNotFoundError(
            f"No default config for {config_name!r} at {default_config} and "
            f"no user override at {user_config}."
        )
    return str(default_config)


def select_env(env_name: str, plugin: str) -> str:
    if USER_ENVS_DIR is not None:
        user_env_dir = Path(USER_ENVS_DIR).expanduser().resolve()
        user_env = user_env_dir / plugin / f"{env_name}.yaml"

        if user_env.exists():
            return str(user_env)

        message = (
            f"Configured user_env_dir is set to {user_env_dir}, but no user "
            f"environment file was found for plugin {plugin!r} and environment "
            f"{env_name!r} at {user_env}. User environment files must be placed "
            "at <user_env_dir>/<plugin>/<env_name>.yaml, must use the .yaml "
            f"extension, and must be named exactly {env_name}.yaml. Expected "
            f"path for this environment: {user_env}. Falling back to the "
            "default environment file."
        )
        smk_logger.warning(message)

    if plugin == "core":
        DEFAULT_ENV_DIR = CORE_ENVS_DIR
    else:
        DEFAULT_ENV_DIR = PLUGINS_DIR / plugin / "envs"

    default_env = DEFAULT_ENV_DIR / f"{env_name}.yaml"
    if not default_env.exists():
        raise FileNotFoundError(
            f"No default env for {env_name!r} at {default_env} and "
            f"no user override at 'user_env_dir'."
        )
    return str(default_env)


rule stage_config:
    input:
        lambda wc: select_config(wc.config_name, wc.plugin),
    output:
        OUTDIRS["configs"] + "/{plugin}/{config_name}.yml",
    run:
        _stage(input[0], output[0], mode="copy")


rule stage_env:
    input:
        lambda wc: select_env(wc.env_name, wc.plugin),
    output:
        OUTDIRS["envs"] + "/{plugin}/{env_name}.yaml",
    run:
        _stage(input[0], output[0], mode="copy")
