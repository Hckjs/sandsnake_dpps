STEREO_COMBINER_DIR = PLUGINS_DIR / "stereo_combiner"
STEREO_COMBINER_SCRIPTS_DIR = STEREO_COMBINER_DIR / "scripts"

test_config_param = Path(config.get(""))

STEREO_COMBINER_OUTDIR = Path(OUTDIRS["plugins"]) / "stero_combiner"
STEREO_COMBINER_PATHS = {
    "merged_source_significances": FERMI_OUTDIR / "merged_source_catalogs.h5",
    "template:catalog_out_dirs": FERMI_OUTDIR / "{catalog}",
    "template:sources": FERMI_OUTDIR / "{catalog}" / "{source}" / "{source}.ecsv",
    "template:source_significances": FERMI_OUTDIR
    / "{catalog}"
    / "{source}"
    / "{source}_significances.ecsv",
}

paths_update("stereo_combiner", STEREO_COMBINER_PATHS)
