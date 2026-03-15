PROFILE ?= ./sandsnake_dpps/profiles/local
CONFIG ?= ./examples/core_analysis_config.yaml
SNAKEFLAGS ?= --printshellcmds --show-failed-logs
BUILD_DIR ?= build


SNAKEFILE = --snakefile sandsnake_dpps/workflow/Snakefile
PROFILEFLAG = --profile $(PROFILE)
CONFIGFLAG = --configfile $(CONFIG)

THESIS_CONFIG_ROOT = examples/thesis

all: | $(BUILD_DIR)
	snakemake $(SNAKEFILE) $(PROFILEFLAG) $(CONFIGFLAG) $(SNAKEFLAGS) --config build_dir=$(BUILD_DIR)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

thesis: | $(BUILD_DIR)
	@set -euo pipefail; \
	echo "=== [1/2] Full Array ==="; \
	snakemake $(SNAKEFILE) $(PROFILE) \
		--config build_dir="$(BUILD_DIR)" \
		--configfile "$(THESIS_CONFIG_ROOT)/full_array/core_analysis_config.yaml" \
		$(SNAKEFLAGS); \
	echo "=== [2/2] Subarrays ==="; \
	for subarray_dir in "$(THESIS_CONFIG_ROOT)"/subarrays/*; do \
		[ -d "$$subarray_dir" ] || continue; \
		echo "--- Subarray: $$subarray_dir ---"; \
		snakemake $(SNAKEFILE) $(PROFILE) \
			--config build_dir="$(BUILD_DIR)" \
			--configfile "$$subarray_dir/core_analysis_config.yaml" \
			$(SNAKEFLAGS); \
	done

# Removes build directory
clean:
	rm -rf $(BUILD_DIR)

.PHONY: all clean
