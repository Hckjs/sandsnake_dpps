# infrastructure
include: "../common/paths.smk"
include: "../common/targets.smk"
include: "../common/providers.smk"
include: "../common/staging.smk"
# core rules
include: "rules/mc/simtel.smk"
include: "rules/mc/dl1.smk"
include: "rules/mc/train.smk"
include: "rules/mc/dl2.smk"
include: "rules/mc/irfs.smk"


def check_target(target):
    # check, dass die benötigten inputs für targets vorhanden sind
    # bzw. die benötigten config einträge (pointing, obstime, etc)
    #
    # Hier checken, ob core target existiert, anstatt im resolver
    #
    # input != target
    pass


# Core target resolver
def resolve_core_targets():
    targets = TARGETS_ENVS()
    config_targets = config.get("targets", {})

    if isinstance(config_targets, list):
        enabled_targets = {target: True for target in config_targets}
    elif isinstance(config_targets, dict):
        enabled_targets = config_targets
    else:
        raise ValueError(
            "config['targets'] must be a mapping of target names to booleans or a list"
        )

    for t, enabled in enabled_targets.items():
        if not enabled:
            continue
        if t == "mc_simtel":
            check_target(t)
            targets.append(TARGETS_MC_SIMTEL("core", resolve=True))
        elif t == "mc_dl1":
            check_target(
                t
            )  # Hier im check möglich, dass input und output mc_dl1 sind -> reprocess
            targets.append(TARGETS_MC_DL1("core", resolve=True))
        elif t == "mc_dl1_merged":
            check_target(t)
            targets.append(TARGETS_MC_DL1_MERGED("core", resolve=True))
        elif t == "rf_energy_regressor":
            check_target(t)
            targets.append(TARGETS_RF_ENERGY_REGRESSOR("core", resolve=True))
        elif t == "rf_particle_classifier":
            check_target(t)
            targets.append(TARGETS_RF_PARTICLE_CLASSIFIER("core", resolve=True))
        elif t == "rf_geometry_reconstructor":
            check_target(t)
            targets.append(TARGETS_RF_GEOMETRY_RECONSTRUCTOR("core", resolve=True))
        elif t == "mc_dl2":
            check_target(t)
            targets.append(TARGETS_MC_DL2("core", resolve=True))
        elif t == "rf_performance_plots":
            check_target(t)
            targets.append(TARGETS_RF_PERFORMANCE_PLOTS("core", resolve=True))
        elif t == "irfs":
            check_target(t)
            targets.append(TARGETS_CUTS("core", resolve=True))
            targets.append(TARGETS_IRFS("core", resolve=True))
            targets.append(TARGETS_BENCHMARKS("core", resolve=True))
        elif t == "irfs_plots":
            check_target(t)
            targets.append(TARGETS_IRF_PLOTS("core", resolve=True))
        else:
            raise ValueError(f"Unknown core target '{t}'")

    return targets


CORE_TARGETS = resolve_core_targets()
