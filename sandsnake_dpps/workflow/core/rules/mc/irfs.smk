
rule optimize_cuts:
    output:
        cuts=OUTPATHS["cuts"],
    input:
        gammas=bind_wildcards(
            mc_dl2_provider, particle="gamma_diffuse", split="test_cuts"
        ),
        protons=bind_wildcards(mc_dl2_provider, particle="proton", split="test_cuts"),
        electrons=bind_wildcards(
            mc_dl2_provider, particle="electron", split="test_cuts"
        ),
        config=PATHS["core:config:optimize_cuts"],
    conda:
        select_env("ctapipe", "core")
    log:
        log=log_path(OUTPATHS["cuts"], ".log"),
        provenance=log_path(OUTPATHS["cuts"], ".provenance"),
    benchmark:
        log_path(OUTPATHS["cuts"], ".benchmark")
    resources:
        mem_mb=4000,
    shell:
        """
        ctapipe-optimize-event-selection \
            --config={input.config} \
            --gamma-file={input.gammas} \
            --proton-file={input.protons} \
            --electron-file={input.electrons} \
            --output={output.cuts} \
            --EventSelectionOptimizer.obs_time="{wildcards.obstime} hour" \
            --log-file={log.log} \
            --provenance-log={log.provenance} \
            --log-level=DEBUG \
        """


rule compute_irfs:
    output:
        irfs=OUTPATHS["irfs"],
        benchmarks=OUTPATHS["benchmarks"],
    input:
        gammas=bind_wildcards(
            mc_dl2_provider, particle="gamma_diffuse", split="test_irfs"
        ),
        protons=bind_wildcards(mc_dl2_provider, particle="proton", split="test_irfs"),
        electrons=bind_wildcards(
            mc_dl2_provider, particle="electron", split="test_irfs"
        ),
        cuts=cuts_provider,
        config=PATHS["core:config:compute_irf"],
    conda:
        select_env("ctapipe", "core")
    log:
        log=log_path(OUTPATHS["irfs"], ".log"),
        provenance=log_path(OUTPATHS["irfs"], ".provenance"),
    benchmark:
        log_path(OUTPATHS["irfs"], ".benchmark")
    resources:
        mem_mb=4000,
    shell:
        """
        ctapipe-compute-irf \
            --config={input.config} \
            --cuts={input.cuts} \
            --gamma-file={input.gammas} \
            --proton-file={input.protons} \
            --electron-file={input.electrons} \
            --output={output.irfs} \
            --benchmark-output={output.benchmarks} \
            --IrfTool.obs_time="{wildcards.obstime} hour" \
            --log-file={log.log} \
            --provenance-log={log.provenance} \
            --log-level=DEBUG \
        """


rule plot_irfs:
    output:
        OUTPATHS["plot_irfs"],
    input:
        irfs=irfs_provider,
        benchmarks=benchmarks_provider,
        cuts=cuts_provider,
        script=CORE_SCRIPTS_DIR / "mc/plot_irfs.py",
        dependency=CORE_SCRIPTS_DIR / "mc/irf_plots.py",
    conda:
        select_env("plotting", "core")
    benchmark:
        log_path(OUTPATHS["plot_irfs"], ".benchmark")
    resources:
        mem_mb=2000,
    shell:
        """
        PYTHONPATH={CORE_SCRIPTS_DIR.parent} \
        python -m scripts.mc.plot_irfs \
        --irfs-file {input.irfs} \
        --cuts-file {input.cuts} \
        --benchmark-file {input.benchmarks} \
        -o {output} \
        """
