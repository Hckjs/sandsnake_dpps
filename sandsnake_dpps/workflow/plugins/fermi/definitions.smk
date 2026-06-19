FERMI_DIR = PLUGINS_DIR / "fermi"
FERMI_SCRIPTS_DIR = FERMI_DIR / "scripts"
REPO_ROOT = WORKFLOW_DIR.parents[1]

DEFAULT_FERMI_CATALOGS_DIR = REPO_ROOT / "resources" / "FERMI_LAT"
DEFAULT_GAMMAPY_DATA_DIR = REPO_ROOT / "resources" / "gammapy-data"


def require_loaded_gammapy_data(gammapy_data_dir: Path):
    """Fail early when the gammapy-data submodule is not checked out."""
    required_paths = [
        gammapy_data_dir / "README.md",
        gammapy_data_dir / "ebl",
    ]
    if all(path.exists() for path in required_paths):
        return

    missing = ", ".join(str(path) for path in required_paths if not path.exists())
    raise ValueError(
        "The Fermi plugin requires the gammapy-data submodule to be loaded at "
        f"{gammapy_data_dir}. Missing: {missing}. Run "
        "`git submodule update --init --recursive resources/gammapy-data` "
        "from the repository root before using the Fermi plugin."
    )


FERMI_CATALOGS_IN_DIR = Path(
    config.get("fermi_catalogs_dir") or DEFAULT_FERMI_CATALOGS_DIR
)
GAMMAPY_DATA_DIR = Path(config.get("gammapy_data_dir") or DEFAULT_GAMMAPY_DATA_DIR)
require_loaded_gammapy_data(GAMMAPY_DATA_DIR)
os.environ["GAMMAPY_DATA"] = str(GAMMAPY_DATA_DIR)

FERMI_CATALOGS = {
    "FGL": FERMI_CATALOGS_IN_DIR / "4FGL_DR4.fit",
    "FHL": FERMI_CATALOGS_IN_DIR / "3FHL.fit",
    "LAC": FERMI_CATALOGS_IN_DIR / "4LAC_DR3_H.fits",
}

FERMI_OUTDIR = OUTDIRS["plugins"] + "/fermi"
FERMI_PATHS = {
    "merged_source_significances": FERMI_OUTDIR + "/merged_source_catalogs.h5",
    "template:catalog_out_dirs": FERMI_OUTDIR + "/{catalog}",
    "template:sources": FERMI_OUTDIR + "/{catalog}/{source}/{source}.ecsv",
    "template:source_significances": FERMI_OUTDIR
    + "/{catalog}/{source}/{source}_significances.ecsv",
}

paths_update("fermi", FERMI_PATHS)

PROCESS_CATALOG_CONFIG = config.get("process_catalog", {})
REDSHIFT_PRIOR_CONFIG = config.get("redshift_priors", {})
CATALOG_CHUNK_CONFIG = config.get("catalog_chunks", {})


def list_source_files(catalog_dir: Path):
    return sorted(
        f
        for f in catalog_dir.glob("*/*.ecsv")
        if f.is_file() and not f.name.endswith("_significances.ecsv")
    )


def get_source_files_from_checkpoint(catalog: str):
    cp = checkpoints.process_catalog.get(catalog=catalog)
    catalog_dir = Path(cp.output.priors).parent
    return list_source_files(catalog_dir)


def source_significance_files_from_sources(source_files):
    return [p.with_name(f"{p.stem}_significances.ecsv") for p in source_files]


def fermi_source_provider(wc):
    fermi_inputs = config.get("fermi_inputs", None).get("processed_sources")
    base = Path(fermi_inputs or FERMI_OUTDIR).expanduser().resolve()
    return str(base / wc.catalog / wc.source / f"{wc.source}.ecsv")


def fermi_priors_provider(wc):
    fermi_inputs = config.get("fermi_inputs", None).get("processed_sources")
    base = Path(fermi_inputs or FERMI_OUTDIR).expanduser().resolve()
    return str(base / wc.catalog / "redshift_priors.ecsv")


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
            Path(FERMI_OUTDIR) / f"{catalog}{str(f).split(f'/{catalog}', 1)[1]}"
            for f in source_files
        ]
        sigs = source_significance_files_from_sources(new_base_source_files)
        return list(map(str, sigs))

    return _call
