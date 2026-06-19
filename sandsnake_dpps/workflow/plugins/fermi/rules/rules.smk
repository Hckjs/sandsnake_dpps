envvars:
    "PYTHONPATH",
    "GAMMAPY_DATA",


checkpoint process_catalog:
    priority: 100
    output:
        priors=PATHS["fermi:template:catalog_out_dirs"] + "/redshift_priors.ecsv",
    input:
        fgl=FERMI_CATALOGS["FGL"],
        lac=FERMI_CATALOGS["LAC"],
        fhl=FERMI_CATALOGS["FHL"],
        script=FERMI_SCRIPTS_DIR / "process_catalog.py",
    params:
        outdir=lambda wc: PATHS["fermi:template:catalog_out_dirs"].format(
            catalog=wc.catalog
        ),
        start=PROCESS_CATALOG_CONFIG.get("start", "2027-04-01"),
        end=PROCESS_CATALOG_CONFIG.get("end", "2028-04-01"),
        step_minutes=PROCESS_CATALOG_CONFIG.get("step_minutes", 20.0),
        alt_min=PROCESS_CATALOG_CONFIG.get("alt_min", 60.0),
        write_plots=PROCESS_CATALOG_CONFIG.get("write_plots", False),
        redshift_lower_quantile=REDSHIFT_PRIOR_CONFIG.get("lower_quantile", 0.16),
        redshift_upper_quantile=REDSHIFT_PRIOR_CONFIG.get("upper_quantile", 0.84),
        redshift_min_sources_per_group=REDSHIFT_PRIOR_CONFIG.get(
            "min_sources_per_group",
            50,
        ),
    wildcard_constraints:
        catalog="4FGL_DR4|3FHL_DR3",
    conda:
        select_env("plotting", "core")
    resources:
        mem_mb=1000,
    script:
        FERMI_SCRIPTS_DIR / "process_catalog.py"


# compare with PROD5 IRFs (LST sub + full array)
rule calc_significances:
    output:
        PATHS["fermi:template:source_significances"],
    input:
        source=fermi_source_provider,
        priors=fermi_priors_provider,
        irfs=TARGETS_IRFS("core", resolve=True),
        benchmarks=TARGETS_BENCHMARKS("core", resolve=True),
    conda:
        select_env("plotting", "core")
    resources:
        mem_mb=100,
    script:
        FERMI_SCRIPTS_DIR / "calc_significances.py"


rule merge_sources:
    output:
        PATHS["fermi:merged_source_significances"],
    input:
        fgl_sources=fermi_source_significance_provider("4FGL_DR4"),
        fhl_sources=fermi_source_significance_provider("3FHL_DR3"),
    conda:
        select_env("plotting", "core")
    resources:
        mem_mb=5000,
    log:
        log_path(PATHS["fermi:merged_source_significances"], ".log"),
    script:
        FERMI_SCRIPTS_DIR / "merge_sources.py"
