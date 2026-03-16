
rule mc_dl1:
    output:
        OUTPATHS["mc_dl1"],
    input:
        data=mc_simtel_provider,
        config=PATHS["core:config:process"],
    conda:
        select_env("ctapipe", "core")
    log:
        provenance=log_path(OUTPATHS["mc_dl1"], ".provenance"),
    resources:
        mem_mb=1500,
    shell:
        """
        ctapipe-process \
        --input {input.data} \
        --output {output} \
        --config {input.config} \
        --provenance-log {log.provenance} \
        --progress \
        """


def symlink_dl1(file_list: list, stamp_path: str) -> list:
    symlink_dir = Path(stamp_path).parent
    symlink_dir.mkdir(parents=True, exist_ok=True)

    for source in file_list:
        src = Path(source)
        if not src.exists():
            print(f"[WARN] source missing: {src}")
            sys.exit(1)

        dst = symlink_dir / src.name
        if dst.is_symlink() or dst.exists():
            dst.unlink()

        rel_src = os.path.relpath(src.resolve(), start=symlink_dir.resolve())
        dst.symlink_to(rel_src)

    Path(stamp_path).touch()


rule mc_dl1_split:
    output:
        OUTPATHS["mc_dl1_split"] + "/.stamp",
    input:
        files=mc_dl1_split_provider,
    resources:
        mem_mb=100,
    run:
        symlink_dl1(input.files, output[0])


rule mc_dl1_merge:
    output:
        OUTPATHS["mc_dl1_merged"],
    input:
        split_stamp=OUTPATHS["mc_dl1_split"] + "/.stamp",
        config=PATHS["core:config:merge"],
    conda:
        select_env("ctapipe", "core")
    params:
        split_dir=OUTPATHS["mc_dl1_split"],
        pattern="*.dl1.h5",
    log:
        log=log_path(OUTPATHS["mc_dl1_merged"], ".log"),
        provenance=log_path(OUTPATHS["mc_dl1_merged"], ".provenance"),
    benchmark:
        log_path(OUTPATHS["mc_dl1_merged"], ".benchmark")
    resources:
        mem_mb=3000,
        slurm_partition="long",
    shell:
        """
        ctapipe-merge \
        --input-dir={params.split_dir}\
        --pattern={params.pattern} \
        --output {output} \
        --log-file {log.log} \
        --provenance-log {log.provenance} \
        --config {input.config} \
        """
