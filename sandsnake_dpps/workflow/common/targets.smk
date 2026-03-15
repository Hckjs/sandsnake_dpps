import os
import glob
import re
import sys
from pathlib import Path

pointings = config.get("pointings")
obstimes = config.get("obstimes")
particles = ["gamma_diffuse", "proton", "electron", "gamma"]
test_splits = ["test_irfs", "test_cuts"]
train_splits = ["train_en", "train_cl_disp"]

valid_particle_splits = [
    # gamma_diffuse
    ("gamma_diffuse", "train_en"),
    ("gamma_diffuse", "train_cl_disp"),
    ("gamma_diffuse", "test_cuts"),
    ("gamma_diffuse", "test_irfs"),
    # proton
    ("proton", "train_cl_disp"),
    ("proton", "test_cuts"),
    ("proton", "test_irfs"),
    # gamma
    ("gamma", "test_cuts"),
    ("gamma", "test_irfs"),
    # electron
    ("electron", "test_cuts"),
    ("electron", "test_irfs"),
]


def get_filenames(zen: int, az: int, particle: str):
    def _collect(base: str, pattern: str):
        files = sorted(glob.glob(os.path.join(base, pattern)))
        return [Path(f).stem.split(".")[0] for f in files]

    mc_dl1 = _ext("mc_dl1")
    mc_simtel = _ext("mc_simtel")

    if mc_dl1:
        base = f"{mc_dl1}/zen_{zen}/az_{az}/{particle}"
        return _collect(base, "*.DL1.h5")

    if mc_simtel:
        base = f"{mc_simtel}/zen_{zen}/az_{az}/{particle}"
        return _collect(base, "*.simtel.gz")

    # TODO: Tristan: Filepattern aus CORSIKA batches/run_ids irgendwie?
    return []


def TARGETS_MC_SIMTEL(base: str, resolve: bool):
    # TODO: Tristan
    target_path = (
        PATHS[f"{base}:mc_simtel"]
        + PATHS["core:template:mc_dl1"]
        + "/{filename}.simtel.gz"
    )
    if resolve:
        return [
            target_path.format(zen=zen, az=az, particle=particle, filename=filename)
            for zen, az in pointings
            for particle in particles
            for filename in get_filenames(zen, az, particle)
        ]
    return target_path


def TARGETS_MC_DL1(base: str, resolve: bool):
    target_path = (
        PATHS[f"{base}:mc_dl1"] + PATHS["core:template:mc_dl1"] + "/{filename}.dl1.h5"
    )
    if resolve:
        return [
            target_path.format(zen=zen, az=az, particle=particle, filename=filename)
            for zen, az in pointings
            for particle in particles
            for filename in get_filenames(zen, az, particle)
        ]
    return target_path


def TARGETS_MC_DL1_MERGED(base: str, resolve: bool):
    target_path = PATHS[f"{base}:mc_dl1_merged"] + PATHS["core:template:mc_dl1_merged"]
    if resolve:
        return [
            target_path.format(zen=zen, az=az, particle=particle, split=split)
            for zen, az in pointings
            for particle, split in valid_particle_splits
        ]

    return target_path


def TARGETS_RF_ENERGY_REGRESSOR(base: str, resolve: bool):
    target_path = PATHS[f"{base}:models"] + PATHS["core:template:rf_energy_regressor"]
    if resolve:
        return [target_path.format(zen=zen, az=az) for zen, az in pointings]

    return target_path


def TARGETS_RF_PARTICLE_CLASSIFIER(base: str, resolve: bool):
    target_path = (
        PATHS[f"{base}:models"] + PATHS["core:template:rf_particle_classifier"]
    )
    if resolve:
        return [target_path.format(zen=zen, az=az) for zen, az in pointings]

    return target_path


def TARGETS_RF_GEOMETRY_RECONSTRUCTOR(base: str, resolve: bool):
    target_path = (
        PATHS[f"{base}:models"] + PATHS["core:template:rf_geometry_reconstructor"]
    )
    if resolve:
        return [target_path.format(zen=zen, az=az) for zen, az in pointings]

    return target_path


def TARGETS_MC_DL2(base: str, resolve: bool):
    target_path = PATHS[f"{base}:mc_dl2"] + PATHS["core:template:mc_dl2"]
    if resolve:
        return [
            target_path.format(zen=zen, az=az, particle=particle, split=split)
            for zen, az in pointings
            for particle in particles
            for split in test_splits
        ]

    return target_path


def TARGETS_RF_PERFORMANCE_PLOTS(base: str, resolve: bool):
    target_path = PATHS[f"{base}:plots"] + PATHS["core:template:rf_performance_plots"]
    if resolve:
        return [target_path.format(zen=zen, az=az) for zen, az in pointings]

    return target_path


def TARGETS_CUTS(base: str, resolve: bool):
    target_path = PATHS[f"{base}:irfs"] + PATHS["core:template:cuts"]
    if resolve:
        return [
            target_path.format(zen=zen, az=az, obstime=obstime)
            for zen, az in pointings
            for obstime in obstimes
        ]

    return target_path


def TARGETS_IRFS(base: str, resolve: bool):
    target_path = PATHS[f"{base}:irfs"] + PATHS["core:template:irfs"]
    if resolve:
        return [
            target_path.format(zen=zen, az=az, obstime=obstime)
            for zen, az in pointings
            for obstime in obstimes
        ]

    return target_path


def TARGETS_BENCHMARKS(base: str, resolve: bool):
    target_path = PATHS[f"{base}:irfs"] + PATHS["core:template:benchmarks"]
    if resolve:
        return [
            target_path.format(zen=zen, az=az, obstime=obstime)
            for zen, az in pointings
            for obstime in obstimes
        ]

    return target_path


def TARGETS_IRF_PLOTS(base: str, resolve: bool):
    target_path = PATHS[f"{base}:plots"] + PATHS["core:template:irf_plots"]
    if resolve:
        return [
            target_path.format(zen=zen, az=az, obstime=obstime)
            for zen, az in pointings
            for obstime in obstimes
        ]

    return target_path


# ------------------------------------------------------------------
# Misc Targets (internal, intermediate, etc.)
# ------------------------------------------------------------------


def TARGETS_ENVS():
    envs = []
    for key, value in PATHS.items():
        parts = key.split(":")
        if key.split(":")[1] == "env":
            envs.append(value)

    return envs
