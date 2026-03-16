import glob
import os
import re
import shutil
import sys
from pathlib import Path
from snakemake.logging import logger as smk_logger


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


def check_targets(enabled_targets):
    inputs = config.get("inputs", {}) or {}
    provided_inputs = {
        key: Path(value).expanduser().resolve()
        for key, value in inputs.items()
        if value not in (None, "", "null")
    }

    if not provided_inputs:
        smk_logger.info(
            "No external inputs configured. New Monte Carlo data will be generated with simtools."
        )
        return

    missing_paths = [
        f"{key} -> {path}" for key, path in provided_inputs.items() if not path.exists()
    ]
    if missing_paths:
        raise ValueError(
            "Configured input path(s) do not exist:\n- " + "\n- ".join(missing_paths)
        )

    active_targets = [target for target, enabled in enabled_targets.items() if enabled]

    conflicting_targets = sorted(
        (set(provided_inputs) & set(active_targets)) - {"mc_dl1"}
    )
    if conflicting_targets:
        raise ValueError(
            "Configured input and output target must not be the same for: "
            + ", ".join(conflicting_targets)
        )

    if "mc_dl1" in provided_inputs and "mc_dl1" in active_targets:
        logger.info(
            "Input and output target are both 'mc_dl1'. DL1 will be reprocessed."
        )

    if any(
        target in {"irfs", "irfs_plots"} for target in active_targets
    ) and not config.get("obstimes"):
        raise ValueError(
            "Configured targets 'irfs'/'irfs_plots' require non-empty 'obstimes'"
        )

    if any(target != "mc_simtel" for target in active_targets) and not config.get(
        "pointings"
    ):
        raise ValueError("Configured targets require non-empty 'pointings'")


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
    if not any(enabled_targets.values()):
        raise ValueError("At least one core target must be set to true")

    check_targets(enabled_targets)

    for t, enabled in enabled_targets.items():
        if not enabled:
            continue
        if t == "mc_simtel":
            targets.append(TARGETS_MC_SIMTEL("core", resolve=True))
        elif t == "mc_dl1":
            targets.append(TARGETS_MC_DL1("core", resolve=True))
        elif t == "mc_dl1_merged":
            targets.append(TARGETS_MC_DL1_MERGED("core", resolve=True))
        elif t == "rf_energy_regressor":
            targets.append(TARGETS_RF_ENERGY_REGRESSOR("core", resolve=True))
        elif t == "rf_particle_classifier":
            targets.append(TARGETS_RF_PARTICLE_CLASSIFIER("core", resolve=True))
        elif t == "rf_geometry_reconstructor":
            targets.append(TARGETS_RF_GEOMETRY_RECONSTRUCTOR("core", resolve=True))
        elif t == "mc_dl2":
            targets.append(TARGETS_MC_DL2("core", resolve=True))
        elif t == "rf_performance_plots":
            targets.append(TARGETS_RF_PERFORMANCE_PLOTS("core", resolve=True))
        elif t == "irfs":
            targets.append(TARGETS_CUTS("core", resolve=True))
            targets.append(TARGETS_IRFS("core", resolve=True))
            targets.append(TARGETS_BENCHMARKS("core", resolve=True))
        elif t == "irfs_plots":
            targets.append(TARGETS_IRF_PLOTS("core", resolve=True))
        else:
            raise ValueError(f"Unknown core target '{t}'")

    return targets


CORE_TARGETS = resolve_core_targets()
