plotting_env = ENVS["plotting"]
ctapipe_env = ENVS["ctapipe"]

scripts_stereo = SCRIPTS["stereo_combiner"]
plots = OUTDIRS["plots"]

stereo_comb = OUTDIRS["stereo_combiner"]
models = OUTDIRS["models"]
mc = OUTDIRS["mc"]
logs = stereo_comb / "logs"
benchmarks = stereo_comb / "benchmarks"


rule stereo_combiner:
    input:
        theta_plot=stereo_comb
        / "plots/irfs/zen_20/az_0/combiner_theta2_reco_lon_lat_zen_20_az_0_ac_full_array.pdf",
        sens_plot=stereo_comb
        / "plots/irfs/zen_20/az_0/combiner_irfs_sens_zen_20_az_0_ac_full_array.pdf",


rule reapply_stereo_combiner:
    output:
        data=stereo_comb
        / "mc/zen_20/az_0/{particle}/{split}/{particle}_zen_20_az_0_{split}_ac_full_array_{combiner}.dl2.h5",
    input:
        data=mc
        / "zen_20/az_0/{particle}/{split}/{particle}_zen_20_az_0_{split}_ac_full_array.dl2.h5",
        script_dependency_1=scripts_stereo / "stereo_combiner.py",
        script_dependency_2=scripts_stereo / "telescope_event_handling.py",
        script=scripts_stereo / "apply_stereo_combiner.py",
    params:
        combiner=lambda wildcards: wildcards.combiner,
    conda:
        ctapipe_env
    wildcard_constraints:
        split="test_irfs|test_cuts|train_cl_disp|train_en",
    benchmark:
        benchmarks / "dl2/{particle}/{split}/benchmark_reapply_stereo_combiner_{particle}_zen_20_az_0_{split}_ac_full_array_{combiner}.txt"
    log:
        logs
        / "dl2/{particle}/{split}/snakemake_reapply_stereo_combiner_{particle}_zen_20_az_0_{split}_ac_full_array_{combiner}.log",
    resources:
        mem_mb=10000,
    shell:
        """
        python -m scripts.stereo_combiner.apply_stereo_combiner \
            --input {input.data}  \
            --output {output.data} \
            --combiner {params.combiner} \
        """


rule stereo_comb_optimize_cuts:
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


rule stereo_comb_create_irfs:
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


rule stereo_comb_plot_theta2_reco_lon_lat:
    output:
        stereo_comb
        / "plots/irfs/zen_20/az_0/combiner_theta2_reco_lon_lat_zen_20_az_0_ac_full_array.pdf",
    input:
        script=scripts_stereo / "plot_stereo_combiner_theta2_reco_lon_lat.py",
        dependency=scripts_stereo / "stereo_combiner_plots.py",
        data=DL2_STEREO_COMB_DATA,
        rc=MATPLOTLIBRC,
    conda:
        plotting_env
    benchmark:
        benchmarks / "plots/benchmark_stereo_comb_plot_theta2_reco_lon_lat_zen_20_az_0_ac_full_array.txt"
    log:
        logs
        / "plots/stereo_comb_plot_theta2_reco_lon_lat_zen_20_az_0_ac_full_array.log",
    resources:
        mem_mb=10000,
    shell:
        """
        python -m scripts.stereo_combiner.plot_stereo_combiner_theta2_reco_lon_lat \
        --input {input.data} \
        --output {output} \
        """


rule stereo_comb_plot_irfs_sens:
    output:
        stereo_comb
        / "plots/irfs/zen_20/az_0/combiner_irfs_sens_zen_20_az_0_ac_full_array.pdf",
    input:
        irfs=IRFS_STEREO_COMB_DATA(file_type="irfs"),
        benchmarks=IRFS_STEREO_COMB_DATA(file_type="benchmarks"),
        script=scripts_stereo / "plot_stereo_combiner_irfs_sens.py",
        dependency=scripts_stereo / "stereo_combiner_plots.py",
        rc=MATPLOTLIBRC,
    conda:
        plotting_env
    benchmark:
        benchmarks / "plots/benchmark_stereo_comb_plot_irfs_sens_zen_20_az_0_ac_full_array.txt"
    log:
        logs / "plots/stereo_comb_plot_irfs_sens_zen_20_az_0_ac_full_array.log",
    resources:
        mem_mb=2000,
    shell:
        """
        python -m scripts.stereo_combiner.plot_stereo_combiner_irfs_sens \
        --input_irfs {input.irfs} \
        --input_benchmarks {input.benchmarks} \
        --output {output} \
        """
