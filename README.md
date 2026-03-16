<p align="center">
  <img src="resources/logo/logo.png" alt="SandSnake DPPS" width="220">
</p>
<h1 align="center">SandSnake DPPS</h1>

## 1. Project Overview

`sandsnake_dpps` is a sandbox Snakemake-based pipeline for processing and analyzing CTAO/ctapipe data, intended as a lightweight playground to test and prototype DPPS workflow components, validate new algorithms, and experiment with new ctapipe functionalities on simulated (TODO:or observational) data without implying a production-ready system. The goal is to provide reproducible analysis workflows for simulated (TODO:and observational) data with configurable input and output targets.

Core building blocks in this repository:

- **Entry point:** `sandsnake_dpps/workflow/Snakefile` (connects core workflow and plugins; defines the `rule all` that aggregates all requested targets)
- **Core workflow:** `sandsnake_dpps/workflow/core/core.smk` (import core rules and common workflow components; defines core targets)
  - `sandsnake_dpps/workflow/core/rules/` (core processing rules, e.g., DL1/DL2 production, IRF calculation, RF training, plotting, benchmarking)
  - `sandsnake_dpps/workflow/core/envs/` (conda environment definitions for core rules)
  - `sandsnake_dpps/workflow/core/configs/` (default core rules configs)
  - `sandsnake_dpps/workflow/core/scripts/` (scripts used by core rules)
- **Shared workflow infrastructure (`common/`):**
  - `sandsnake_dpps/workflow/common/paths.smk` (canonical directory and path templates)
  - `sandsnake_dpps/workflow/common/providers.smk` (input providers that choose between `core` outputs and staged external inputs)
  - `sandsnake_dpps/workflow/common/staging.smk` (staging rules that normalize external data into the internal layout)
  - `sandsnake_dpps/workflow/common/targets.smk` (declares the requested end products; Snakemake resolves all required upstream steps automatically)
- **Profiles:** `sandsnake_dpps/profiles/local`, `sandsnake_dpps/profiles/slurm` (default Snakemake profiles for local and SLURM execution)
- **Plugins:** `sandsnake_dpps/plugins/` (example plugin workflows that can be optionally included, e.g., Fermi analysis)
- **Example configurations:**
  - `examples/core_analysis_config.yaml`
  - `examples/fermi_analysis_config.yaml`

---

## 2. How the package works (pipeline logic)

The pipeline is controlled by a YAML configuration and then resolved by Snakemake into a DAG (dependency graph).

1. **Load configuration**
   - In normal usage, you run the workflow through the `Makefile` and set the config via `CONFIG=...` (this is forwarded to Snakemake as `--configfile <file.yaml>`). Default is set to `examples/core_analysis_config.yaml`.
   - Most important fields are:
     - `targets` (which final artifacts should be built)
     - `inputs.*` (source paths for existing data)
     - `plugins` (which optional plugins to include)

2. **Select core targets**
   - Desired core targets are listed in `targets`, for example:
     - `mc_dl1`
     - `rf_performance_plots`
     - `irfs`
     - `irfs_plots`
   - Only requested targets are resolved; unknown names raise an error (`Unknown core target`).

3. **Optionally load plugins**
   - Additional workflows can be included via `plugins`, for example:
     - `fermi`
   - Plugins register their own final artifacts, which are built together with core targets.

4. **Automatically resolve dependencies**
   - Snakemake builds only the artifacts required to satisfy selected targets.

**Understand the `common/` principle (providers + staging)**
   - The workflow uses a two-layer strategy so rules can stay stable even when input sources vary.
   - `providers.smk` acts as a routing layer: if an `inputs.*` path is configured, providers resolve to `stage:*` locations; otherwise they resolve to `core:*` locations generated inside the pipeline.
   - `staging.smk` imports external artifacts (e.g., DL1/DL2, RF models, IRFs) into a normalized internal structure (typically under `build/<run_tag>/staging/...`) so downstream rules can use one consistent path scheme.
   - `paths.smk` centralizes path templates and directory namespaces (`core`, `stage`, `plugin`) to keep file naming and folder layout deterministic.

`rule all` in `sandsnake_dpps/workflow/Snakefile` aggregates both **core targets** and **plugin targets** as joint final outputs.

---

## 3. Prerequisites

- **Conda or Mamba** (environment management)

Available environment definitions:

- `environment.yml` (project-level environment)
- `sandsnake_dpps/workflow/core/envs/*.yaml` (rule-specific default core workflow environments)

### Expected input directory layout and naming

If you provide external data via `inputs.*`, files should follow the workflow's template-based layout under the configured input root.

- **MC DL1 input root (`inputs.mc_dl1`)**
  - Expected directory pattern: `zen_<ZEN>/az_<AZ>/<particle>/`
  - Expected file naming (external): `<filename>.DL1.h5` (uppercase `DL1`), e.g. `run123.DL1.h5`
  - During staging, files are normalized to internal lowercase naming (`.dl1.h5`). (TODO)

- **MC simtel input root (`inputs.mc_simtel`)**
  - TODO

- **Template-driven products (when provided externally)**
  - `inputs.mc_dl1_merged`: `zen_<ZEN>/az_<AZ>/<particle>/<split>/<particle>_zen_<ZEN>_az_<AZ>_<split>_merged.dl1.h5`
  - `inputs.mc_dl2`: `zen_<ZEN>/az_<AZ>/<particle>/<split>/<particle>_zen_<ZEN>_az_<AZ>_<split>.dl2.h5`
  - `inputs.models`: `zen_<ZEN>/az_<AZ>/rf_energy_regressor_zen_<ZEN>_az_<AZ>.pkl`, `rf_particle_classifier_...`, `rf_geometry_reconstructor_...`
  - `inputs.cuts`: `zen_<ZEN>/az_<AZ>/cuts_zen_<ZEN>_az_<AZ>_obs_<OBSTIME>_hours.fits`
  - `inputs.irfs`: `zen_<ZEN>/az_<AZ>/irfs_zen_<ZEN>_az_<AZ>_obs_<OBSTIME>_hours.fits.gz`
  - `inputs.benchmarks`: `zen_<ZEN>/az_<AZ>/benchmarks_zen_<ZEN>_az_<AZ>_obs_<OBSTIME>_hours.fits.gz`

- **Placeholder constraints**
  - `<particle>` must be one of: `gamma_diffuse`, `proton`, `electron`, `gamma`.
  - `<split>` is from: `train_en`, `train_cl_disp`, `test_cuts`, `test_irfs` (valid combinations are particle-dependent in the workflow logic).
  - `<OBSTIME>` must match one of the values listed in config `obstimes` (typically hour values such as `0.5`, `5` or `50`, reflected in filenames as `obs_<OBSTIME>_hours`).
  - `<ZEN>` and `<AZ>` must be integers, matching the values defined in config `pointings` (e.g. `[20, 0]`).


If these naming/layout conventions are not met, staging/provider resolution may succeed only partially or produce missing-target behavior.

### Notes on key parameters in `examples/core_analysis_config.yaml`

The default core config is a good reference for which parameters control behavior:

- `run_tag`: names the run subdirectory under `build/` (e.g., `build/<run_tag>/...`).
- `user_config_dir` and `user_env_dir`: optional override roots for staged config/env files; files must match the existing naming scheme. Core config overrides must be in `core/` and use the same config names as defaults (plugin configs can be placed under `<plugin>/` with matching names). Env override files must also use matching env names (core under `core/`, plugin envs under `<plugin>/`).
- `targets`: defines which core end products should be built (e.g., `mc_dl1`, `irfs`, `irfs_plots`).
- `inputs.*`: external roots for existing artifacts (`mc_dl1`, `mc_dl2`, `models`, `irfs`, etc.); `null` means the pipeline builds that stage internally.
- `pointings`: list of `[ZEN, AZ]` integer pairs used in template expansion and target generation.
- `obstimes`: observation times (hours) used in cut/IRF/benchmark outputs and filenames.
- `train_size` and `cuts_size`: split fractions used by the core rules for training/testing/cut optimization subsets.
- `plugins`: optional plugin list, usually `null` for pure core runs.

---

## 4. Quickstart example: generate sensitivity-related outputs from existing DL1 files (local)

If you want to generate sensitivity-related artifacts/plots from existing DL1 files:

1. Create your own core config based on `examples/core_analysis_config.yaml`.
2. Keep the structure, but set `inputs.mc_dl1` to your local DL1 root path.
3. Ensure your external DL1 files follow the naming and directory scheme from the previous section (**Expected input directory layout and naming**).
4. Run with your config, for example:

Dry run: only builds the workflow DAG and shows which files/rules would be executed (no jobs are actually run, no outputs are created).
```bash
make CONFIG=path/to/your/core_analysis_config.yaml SNAKEFLAGS="--dry-run --printshellcmds --show-failed-logs"
```

Use a custom Snakemake profile (e.g. local/SLURM) by pointing PROFILE to the profile directory; the workflow is executed with the settings defined in that profile.
```bash
make CONFIG=path/to/your/core_analysis_config.yaml PROFILE=path/to/your/profile/dir
```

Important variables from the `Makefile`:

- `CONFIG` – path to the configuration file (default: `examples/core_analysis_config.yaml`)
- `PROFILE` – path to profile dir including a `config.yaml` (default: `sandsnake_dpps/profiles/local`)
- `BUILD_DIR` – output directory (default: `build`)
- `SNAKEFLAGS` – additional Snakemake flags e.g. `--dry-run` (default includes `--printshellcmds --show-failed-logs`)


**Note:** If `mc_dl1` is set both as an external input (`inputs.mc_dl1`) and as a requested target (`targets` contains `mc_dl1`), the workflow reprocesses DL1 (instead of only consuming the staged input).

---

## 5. Troubleshooting

- **`Unknown core target`**
  - Verify target names in `targets` (spelling must match a known core target).

- **Empty or missing outputs**
  - Check all paths under `inputs.*`. Ensure they point to existing data with the expected directory structure and naming conventions.
  - Ensure selected `targets` have all required input data.

- **Conda environment issues**
  - Verify `use-conda` is enabled in the selected profile.
  - Ensure env files are available and correctly named (core envs under `core/envs/`, plugin envs under `<plugin>/envs/`).
  - By default, Snakemake usually creates Conda envs under the repo-local `.snakemake/` directory. If you hit home quota limits on a cluster, set `conda-prefix` in the selected profile config to a different filesystem path with enough space.

- **Path resolution / relative vs absolute paths**
  - Relative paths are generally resolved by the workflow, but to avoid edge cases it is still safer to use absolute paths wherever possible (especially for `inputs.*`, user config/env dirs, and custom output roots).

- **Default core configs are not being overridden**
  - Check that your `user_config_dir` path (from the core config) contains the expected naming scheme.
  - Core config overrides must be located under `<user_config_dir>/core/` and use the same config filenames as the defaults (plugin overrides under `<user_config_dir>/<plugin>/` with matching names).

---

## Command reference

```bash
# Using tmux
tmux new -s ssdpps 'make CONFIG=path/to/your/core_analysis_config.yaml PROFILE=path/to/your/profile/dir | tee -a ~/logs/ssdpps.log'
# detach with Ctrl+B, D
```
```bash
# attach back with
tmux attach -t ssdpps
```
