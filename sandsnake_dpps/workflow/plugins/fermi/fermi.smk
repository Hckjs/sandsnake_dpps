# Fermi plugin entrypoint
if USER_CONFIGS_DIR is None:
    raise ValueError(
        "Missing user_config_dir in core_analysis_config.yaml. "
        "Provide user_config_dir which includes a fermi_config.yaml file."
    )


configfile: Path(USER_CONFIGS_DIR) / "fermi_analysis_config.yaml"


include: "definitions.smk"
include: "rules/rules.smk"


def resolve_fermi_targets():
    fermi_targets = []
    for t in config.get("fermi_targets", []):
        if t == "merged_source_significances":
            fermi_targets.append(PATHS["fermi:merged_source_significances"])
        if t == "processed_sources":
            catalogs = ["4FGL_DR4", "3FHL_DR3"]
            tar = expand(
                FERMI_PATHS["template:catalog_out_dirs"] / ".stamp", catalog=catalogs
            )
            fermi_targets.extend(tar)

    return fermi_targets


FERMI_TARGETS = resolve_fermi_targets()
register_plugin_targets(FERMI_TARGETS)
