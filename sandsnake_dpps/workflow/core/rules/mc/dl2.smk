
rule mc_dl2:
    output:
        OUTPATHS["mc_dl2"],
    input:
        data=mc_dl1_merged_provider,
        model_e_reg=rf_energy_regressor_provider,
        model_p_clf=rf_particle_classifier_provider,
        model_disp=rf_geometry_reconstructor_provider,
        config=PATHS["core:config:apply_models"],
    conda:
        select_env("ctapipe", "core")
    log:
        log=log_path(OUTPATHS["mc_dl2"], ".log"),
        provenance=log_path(OUTPATHS["mc_dl2"], ".provenance"),
    benchmark:
        log_path(OUTPATHS["mc_dl2"], ".benchmark")
    resources:
        slurm_partition="long",
        mem_mb=20000,
    shell:
        """
        ctapipe-apply-models \
            --input {input.data}  \
            --output {output} \
            --log-file {log.log} \
            --provenance-log {log.provenance} \
            --config {input.config} \
            --reconstructor {input.model_e_reg} \
            --reconstructor {input.model_p_clf} \
            --reconstructor {input.model_disp} \
        """


rule plot_rf_performance:
    output:
        OUTPATHS["rf_performance_plots"],
    input:
        model_e_reg=rf_energy_regressor_provider,
        model_p_clf=rf_particle_classifier_provider,
        model_disp=rf_geometry_reconstructor_provider,
        gamma=bind_wildcards(
            mc_dl2_provider, particle="gamma_diffuse", split="test_irfs"
        ),
        proton=bind_wildcards(mc_dl2_provider, particle="proton", split="test_irfs"),
        script=CORE_SCRIPTS_DIR / "mc/plot_dl2_rf_performance.py",
        dependency=CORE_SCRIPTS_DIR / "mc/dl2_rf_performance_plots.py",
        config=PATHS["core:config:train_models"],
    conda:
        select_env("plotting", "core")
    benchmark:
        log_path(OUTPATHS["rf_performance_plots"], ".benchmark")
    resources:
        mem_mb=25000,
    shell:
        """
        python {input.script} \
        --gammas {input.gamma} \
        --protons {input.proton} \
        --disp-model {input.model_disp} \
        --e-reg-model {input.model_e_reg} \
        --clf-model {input.model_p_clf} \
        --config {input.config} \
        --output {output} \
        """
