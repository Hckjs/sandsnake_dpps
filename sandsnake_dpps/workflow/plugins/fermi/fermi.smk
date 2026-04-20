from pathlib import Path
import re

# Fermi plugin entrypoint
if USER_CONFIGS_DIR is None:
    raise ValueError(
        "Missing user_config_dir in core_analysis_config.yaml. "
        "Provide user_config_dir which includes a fermi_analysis_config.yaml file."
    )


configfile: Path(USER_CONFIGS_DIR) / "plugins/fermi/fermi_analysis_config.yaml"


include: "definitions.smk"
include: "rules/rules.smk"


def check_fermi_targets(enabled_targets):
    return


def resolve_fermi_targets():
    fermi_targets = []
    config_targets = config.get("fermi_targets", {})

    if isinstance(config_targets, list):
        enabled_targets = {target: True for target in config_targets}
    elif isinstance(config_targets, dict):
        enabled_targets = config_targets
    else:
        raise ValueError(
            "config['targets'] must be a mapping of target names to booleans or a list"
        )
    if not any(enabled_targets.values()):
        raise ValueError("At least one fermi target must be set to true")

    check_fermi_targets(enabled_targets)

    for t, enabled in enabled_targets.items():
        if not enabled:
            continue
        if t == "merged_source_significances":
            fermi_targets.append(PATHS["fermi:merged_source_significances"])
        if t == "processed_sources":
            catalogs = ["4FGL_DR4", "3FHL_DR3"]
            tar = expand(
                PATHS["fermi:template:catalog_out_dirs"] / ".stamp", catalog=catalogs
            )
            fermi_targets.extend(tar)

    return fermi_targets


FERMI_TARGETS = resolve_fermi_targets()
register_plugin_targets(FERMI_TARGETS)
