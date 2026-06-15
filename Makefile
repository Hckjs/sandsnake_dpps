PROFILE ?= ./sandsnake_dpps/profiles/vollmond
CONFIG ?= ./examples/core_analysis_config.yaml
BUILD_DIR ?= build


SNAKEFILE = --snakefile sandsnake_dpps/workflow/Snakefile
PROFILEFLAG = --profile $(PROFILE)
CONFIGFLAG = --configfile $(CONFIG)
SNAKEFLAGS ?= #--dry-run

THESIS_CONFIG_ROOT = ./examples/thesis

all: | $(BUILD_DIR)
	snakemake $(SNAKEFILE) $(PROFILEFLAG) $(CONFIGFLAG) $(SNAKEFLAGS) \
		--config build_dir=$(BUILD_DIR)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

thesis: | $(BUILD_DIR)
	echo "=== [1/2] Full Array ==="; \
	snakemake $(SNAKEFILE) $(PROFILEFLAG) \
		--config build_dir="$(BUILD_DIR)" \
		--configfile "$(THESIS_CONFIG_ROOT)/configs/full_array/core_analysis_config.yaml" \
		$(SNAKEFLAGS); \
	echo "=== [2/2] Subarrays ==="; \
	for subarray_dir in "$(THESIS_CONFIG_ROOT)"/configs/subarrays/*; do \
		[ -d "$$subarray_dir" ] || continue; \
		echo "--- Subarray: $$subarray_dir ---"; \
		snakemake $(SNAKEFILE) $(PROFILEFLAG) \
			--config build_dir="$(BUILD_DIR)" \
			--configfile "$$subarray_dir/core_analysis_config.yaml" \
			$(SNAKEFLAGS); \
	done

# Removes build directory
clean:
	rm -rf $(BUILD_DIR)

.PHONY: all clean
