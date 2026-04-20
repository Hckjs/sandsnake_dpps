# Stereo Combiner plugin entrypoint
if USER_CONFIGS_DIR is None:
    raise ValueError(
        "Missing user_config_dir in core_analysis_config.yaml. "
        "Provide user_config_dir which includes a stereo_combiner_config.yaml file."
    )


configfile: Path(USER_CONFIGS_DIR) / "plugins/stereo_combiner/stereo_combiner_config.yaml"


include: "definitions.smk"
include: "rules/rules.smk"


def check_stereo_combiner_targets(enabled_targets):
    return


def resolve_stereo_combiner_targets():
    stereo_combiner_targets = []
    config_targets = config.get("stereo_combiner_targets", {})

    if isinstance(config_targets, list):
        enabled_targets = {target: True for target in config_targets}
    elif isinstance(config_targets, dict):
        enabled_targets = config_targets
    else:
        raise ValueError(
            "config['targets'] must be a mapping of target names to booleans or a list"
        )
    if not any(enabled_targets.values()):
        raise ValueError("At least one stereo combiner target must be set to true")

    check_stereo_combiner_targets(enabled_targets)

    for t, enabled in enabled_targets.items():
        if not enabled:
            continue
        if t == "plot_reco_lon_lat":
            stereo_combiner_targets.append(PATHS["stereo_combiner:plot_reco_lon_lat"])
        if t == "plot_irfs":
            stereo_combiner_targets.append(PATHS["stereo_combiner:plot_irfs"])

    return stereo_combiner_targets


STEREO_COMBINER_TARGETS = resolve_stereo_combiner_targets()
register_plugin_targets(STEREO_COMBINER_TARGETS)
