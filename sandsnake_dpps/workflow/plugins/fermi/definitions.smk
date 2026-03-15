from pathlib import Path
import re


FERMI_DIR = PLUGINS_DIR / "fermi"
FERMI_SCRIPTS_DIR = FERMI_DIR / "scripts"

FERMI_CATALOGS_IN_DIR = Path(config.get("fermi_catalogs_dir"))
GAMMAPY_DATA_DIR = Path(config.get("gammapy_data_dir"))
os.environ["GAMMAPY_DATA"] = str(GAMMAPY_DATA_DIR)

FERMI_CATALOGS = {
    "FGL": FERMI_CATALOGS_IN_DIR / "4FGL_DR4.fit",
    "FHL": FERMI_CATALOGS_IN_DIR / "3FHL_DR3.fit",
    "LAC": FERMI_CATALOGS_IN_DIR / "4LAC_DR3.fits",
}

FERMI_OUTDIR = Path(OUTDIRS["plugins"]) / "fermi"
FERMI_PATHS = {
    "merged_source_significances": FERMI_OUTDIR / "merged_source_catalogs.h5",
    "template:catalog_out_dirs": FERMI_OUTDIR / "{catalog}",
    "template:sources": FERMI_OUTDIR / "{catalog}" / "{source}" / "{source}.ecsv",
    "template:source_significances": FERMI_OUTDIR
    / "{catalog}"
    / "{source}"
    / "{source}_significances.ecsv",
}

paths_update("fermi", FERMI_PATHS)


def list_source_files(catalog_dir: Path):
    return sorted(
        f
        for f in catalog_dir.glob("*/*.ecsv")
        if f.is_file() and not f.name.endswith("_significances.ecsv")
    )


def get_source_files_from_checkpoint(catalog: str):
    cp = checkpoints.process_catalog.get(catalog=catalog)
    catalog_dir = Path(cp.output[0])
    return list_source_files(catalog_dir)


def source_significance_files_from_sources(source_files):
    return [p.with_name(f"{p.stem}_significances.ecsv") for p in source_files]


def fermi_source_provider(wc):
    fermi_inputs = config.get("fermi_inputs", None).get("processed_sources")
    base = Path(fermi_inputs or FERMI_OUTDIR).expanduser().resolve()
    return str(base / wc.catalog / wc.source / f"{wc.source}.ecsv")


def fermi_source_significance_provider(catalog: str):
    def _call(wc):
        fermi_inputs = config.get("fermi_inputs", None).get("processed_sources")
        base = Path(fermi_inputs or FERMI_OUTDIR).expanduser().resolve()

        if fermi_inputs:
            source_files = sorted(
                f
                for f in base.glob(f"{catalog}/*/*.ecsv")
                if f.is_file() and not f.name.endswith("_significances.ecsv")
            )
        else:
            source_files = get_source_files_from_checkpoint(catalog)

        new_base_source_files = [
            FERMI_OUTDIR / f"{catalog}{str(f).split(f'/{catalog}', 1)[1]}"
            for f in source_files
        ]
        sigs = source_significance_files_from_sources(new_base_source_files)
        return list(map(str, sigs))

    return _call
