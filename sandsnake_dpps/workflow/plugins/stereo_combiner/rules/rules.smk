rule reapply_stereo_combiner:
    output:
        PATHS["stereo_combiner:mc_dl2"],
    input:
        data=bind_wildcards(mc_dl2_provider, zen="20", az="180"),
        script=STEREO_COMBINER_SCRIPTS_DIR / "apply_stereo_combiner.py",
        sd_1=STEREO_COMBINER_SCRIPTS_DIR / "stereo_combiner.py",
        sd_2=STEREO_COMBINER_SCRIPTS_DIR / "telescope_event_handling.py",
    params:
        combiner=lambda wildcards: wildcards.combiner,
    conda:
        select_env("ctapipe", "core")
    log:
        log=log_path(PATHS["stereo_combiner:mc_dl2"], ".log"),
    benchmark:
        log_path(PATHS["stereo_combiner:mc_dl2"], ".benchmark")
    resources:
        mem_mb=10000,
    shell:
        """
        PYTHONPATH={STEREO_COMBINER_DIR} \
        python -m scripts.apply_stereo_combiner \
            --input {input.data}  \
            --output {output} \
            --combiner {params.combiner} \
        """


rule stereo_comb_optimize_cuts:
    output:
        cuts=PATHS["stereo_combiner:cuts"],
    input:
        gammas=bind_wildcards(
            PATHS["stereo_combiner:mc_dl2"], particle="gamma", split="test_cuts"
        ),
        protons=bind_wildcards(
            PATHS["stereo_combiner:mc_dl2"], particle="proton", split="test_cuts"
        ),
        electrons=bind_wildcards(
            PATHS["stereo_combiner:mc_dl2"], particle="electron", split="test_cuts"
        ),
        config=PATHS["core:config:optimize_cuts"],
    conda:
        select_env("ctapipe", "core")
    log:
        log=log_path(PATHS["stereo_combiner:cuts"], ".log"),
        provenance=log_path(PATHS["stereo_combiner:cuts"], ".provenance"),
    benchmark:
        log_path(PATHS["stereo_combiner:cuts"], ".benchmark")
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
            --log-file={log.log} \
            --provenance-log={log.provenance} \
            --log-level=DEBUG \
        """


rule stereo_comb_create_irfs:
    output:
        irfs=PATHS["stereo_combiner:irfs"],
        benchmarks=PATHS["stereo_combiner:benchmarks"],
    input:
        gammas=bind_wildcards(
            PATHS["stereo_combiner:mc_dl2"], particle="gamma", split="test_irfs"
        ),
        protons=bind_wildcards(
            PATHS["stereo_combiner:mc_dl2"], particle="proton", split="test_irfs"
        ),
        electrons=bind_wildcards(
            PATHS["stereo_combiner:mc_dl2"], particle="electron", split="test_irfs"
        ),
        cuts=PATHS["stereo_combiner:cuts"],
        config=PATHS["core:config:compute_irf"],
    conda:
        select_env("ctapipe", "core")
    log:
        log=log_path(PATHS["stereo_combiner:irfs"], ".log"),
        provenance=log_path(PATHS["stereo_combiner:irfs"], ".provenance"),
    benchmark:
        log_path(PATHS["stereo_combiner:irfs"], ".benchmark")
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
            --log-file={log.log} \
            --provenance-log={log.provenance} \
            --log-level=DEBUG \
        """


rule stereo_comb_plot_theta2_reco_lon_lat:
    output:
        PATHS["stereo_combiner:plot_reco_lon_lat"],
    input:
        data=STEREO_COMBINER_MC_DL2,
        script=STEREO_COMBINER_SCRIPTS_DIR
        / "plot_stereo_combiner_theta2_reco_lon_lat.py",
        sd_1=STEREO_COMBINER_SCRIPTS_DIR / "stereo_combiner_plots.py",
    conda:
        select_env("plotting", "core")
    log:
        log_path(PATHS["stereo_combiner:plot_reco_lon_lat"], ".log"),
    benchmark:
        log_path(PATHS["stereo_combiner:plot_reco_lon_lat"], ".benchmark")
    resources:
        mem_mb=10000,
    shell:
        """
        PYTHONPATH={WORKFLOW_DIR} \
        python -m plugins.stereo_combiner.scripts.plot_stereo_combiner_theta2_reco_lon_lat \
        --input {input.data} \
        --output {output} \
        """


rule stereo_comb_plot_irfs_sens:
    output:
        PATHS["stereo_combiner:plot_irfs"],
    input:
        irfs=STEREO_COMBINER_IRFS(file_type="irfs"),
        benchmarks=STEREO_COMBINER_IRFS(file_type="benchmarks"),
        script=STEREO_COMBINER_SCRIPTS_DIR / "plot_stereo_combiner_irfs_sens.py",
        sd_1=STEREO_COMBINER_SCRIPTS_DIR / "stereo_combiner_plots.py",
    conda:
        select_env("plotting", "core")
    log:
        log_path(PATHS["stereo_combiner:plot_irfs"], ".log"),
    benchmark:
        log_path(PATHS["stereo_combiner:plot_irfs"], ".benchmark")
    resources:
        mem_mb=2000,
    shell:
        """
        PYTHONPATH={STEREO_COMBINER_DIR} \
        python -m scripts.plot_stereo_combiner_irfs_sens \
        --input_irfs {input.irfs} \
        --input_benchmarks {input.benchmarks} \
        --output {output} \
        """
