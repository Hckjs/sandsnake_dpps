
rule train_energy_regressor:
    output:
        OUTPATHS["train_energy_regressor"],
    input:
        gammas=bind_wildcards(
            mc_dl1_merged_provider, particle="gamma_diffuse", split="train_en"
        ),
        config=PATHS["core:config:train_models"],
    conda:
        select_env("ctapipe", "core")
    log:
        log=log_path(OUTPATHS["train_energy_regressor"], ".log"),
        provenance=log_path(OUTPATHS["train_energy_regressor"], ".provenance"),
    benchmark:
        log_path(OUTPATHS["train_energy_regressor"], ".benchmark")
    resources:
        mem_mb=24000,
        slurm_partition="long",
    threads: 8
    shell:
        """
        ctapipe-train-energy-regressor \
            --config {input.config} \
            --input {input.gammas} \
            --output {output} \
            --n-jobs {threads} \
            --log-file {log.log} \
            --provenance-log {log.provenance} \
            --log-level=INFO \
            """


rule apply_energy_regressor_train:
    output:
        OUTPATHS["apply_energy_regressor_train"],
    input:
        data=bind_wildcards(mc_dl1_merged_provider, split="train_cl_disp"),
        model_e_reg=rf_energy_regressor_provider,
        config=PATHS["core:config:apply_models"],
    conda:
        select_env("ctapipe", "core")
    log:
        log=log_path(OUTPATHS["apply_energy_regressor_train"], ".log"),
        provenance=log_path(OUTPATHS["apply_energy_regressor_train"], ".provenance"),
    benchmark:
        log_path(OUTPATHS["apply_energy_regressor_train"], ".benchmark")
    resources:
        mem_mb=15000,
        slurm_partition="long",
    shell:
        """
        ctapipe-apply-models \
            --input {input.data}  \
            --output {output} \
            --log-file {log.log} \
            --provenance-log {log.provenance} \
            --config {input.config} \
            --reconstructor {input.model_e_reg} \
        """


rule train_disp_reconstructor:
    output:
        OUTPATHS["train_disp_reconstructor"],
    input:
        gammas=bind_wildcards(
            OUTPATHS["apply_energy_regressor_train"], particle="gamma_diffuse"
        ),
        config=PATHS["core:config:train_models"],
    conda:
        select_env("ctapipe", "core")
    log:
        log=log_path(OUTPATHS["train_disp_reconstructor"], ".log"),
        provenance=log_path(OUTPATHS["train_disp_reconstructor"], ".provenance"),
    benchmark:
        log_path(OUTPATHS["train_disp_reconstructor"], ".benchmark")
    resources:
        mem_mb=30000,
        slurm_partition="long",
    threads: 8
    shell:
        """
        ctapipe-train-disp-reconstructor \
        --config {input.config} \
        --input {input.gammas} \
        --output {output} \
        --n-jobs {threads} \
        --log-file {log.log} \
        --provenance-log {log.provenance} \
        """


rule train_particle_classifier:
    output:
        OUTPATHS["train_particle_classifier"],
    input:
        gammas=bind_wildcards(
            OUTPATHS["apply_energy_regressor_train"], particle="gamma_diffuse"
        ),
        protons=bind_wildcards(
            OUTPATHS["apply_energy_regressor_train"], particle="proton"
        ),
        config=PATHS["core:config:train_models"],
    conda:
        select_env("ctapipe", "core")
    log:
        log=log_path(OUTPATHS["train_particle_classifier"], ".log"),
        provenance=log_path(OUTPATHS["train_particle_classifier"], ".provenance"),
    benchmark:
        log_path(OUTPATHS["train_particle_classifier"], ".benchmark")
    resources:
        mem_mb=40000,
        slurm_partition="long",
    threads: 8
    shell:
        """
        ctapipe-train-particle-classifier \
        --config {input.config} \
        --signal {input.gammas} \
        --background {input.protons} \
        --output {output} \
        --n-jobs {threads} \
        --log-file {log.log} \
        --provenance-log {log.provenance} \
        """
