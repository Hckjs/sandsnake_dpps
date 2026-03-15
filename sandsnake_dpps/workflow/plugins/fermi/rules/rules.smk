envvars:
    "PYTHONPATH",
    "GAMMAPY_DATA",


ruleorder: reprocess_dl1_test_arrays > stage_mc_dl1_merged


rule reprocess_dl1_test_arrays:
    output:
        PATHS["stage:mc_dl1_merged"] + PATHS["core:template:mc_dl1_merged"],
    input:
        data=_ext("mc_dl1_merged", re_none=False) + PATHS["core:template:mc_dl1_merged"],
        config=PATHS["core:config:process"],
    conda:
        select_env("ctapipe", "core")
    log:
        log=log_path(OUTPATHS["mc_dl1_merged"], ".log"),
        provenance=log_path(OUTPATHS["mc_dl1_merged"], ".provenance"),
    benchmark:
        log_path(OUTPATHS["mc_dl1_merged"], ".benchmark")
    resources:
        mem_mb=24000,
    shell:
        """
        ctapipe-process \
        --input {input.data} \
        --output {output[0]} \
        --log-file {log.log} \
        --provenance-log {log.provenance} \
        --config {input.config} {input.config}\
        --progress \
        """


# TODO: 4FHL recently announced but not yet published.
checkpoint process_catalog:
    output:
        outdir=directory(PATHS["fermi:template:catalog_out_dirs"]),
        stamp=PATHS["fermi:template:catalog_out_dirs"] / ".stamp",
    input:
        fgl=FERMI_CATALOGS["FGL"],
        lac=FERMI_CATALOGS["LAC"],
        fhl=FERMI_CATALOGS["FHL"],
        script=FERMI_SCRIPTS_DIR / "process_catalog.py",
    wildcard_constraints:
        catalog="4FGL_DR4|3FHL_DR3",
    conda:
        select_env("plotting", "core")
    resources:
        mem_mb=1000,
    shell:
        r"""
        set -euo pipefail
        python {input.script} \
            --catalog {wildcards.catalog} \
            --fgl {input.fgl} \
            --lac {input.lac} \
            --fhl {input.fhl} \
            --outdir {output.outdir}
        touch {output.stamp}
        """


# compare with PROD5 IRFs (LST sub + full array)
rule calc_significances:
    output:
        PATHS["fermi:template:source_significances"],
    input:
        source=fermi_source_provider,
        irfs=TARGETS_IRFS("core", resolve=True),
        benchmarks=TARGETS_BENCHMARKS("core", resolve=True),
        stamp=PATHS["fermi:template:catalog_out_dirs"] / ".stamp",
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
