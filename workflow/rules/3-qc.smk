# workflow/rules/3-qc.smk

rule fastp_qc:
    input:
        r1 = rules.extract_chlamydiales_reads.output.r1,
        r2 = rules.extract_chlamydiales_reads.output.r2,
    output:
        html   = OUTDIR / "{sample}" / "qc" / "{sample}_scrub.html",
        json   = OUTDIR / "{sample}" / "qc" / "{sample}_scrub.json",
        status = OUTDIR / "{sample}" / "status" / "fastp_qc.{sample}.txt",
    threads: config["threads"]["fastp"]
    # if fastp is already in trim.yaml, keep this; otherwise point to a qc.yaml with fastp+multiqc
    conda: "../envs/trim.yaml"
    log: OUTDIR / "{sample}" / "log" / "trim.scrub.{sample}.log"
    shell: r"""
      mkdir -p "$(dirname {output.html})" "$(dirname {output.status})" "$(dirname {log})"
      fastp \
        -i {input.r1} \
        -I {input.r2} \
        -h {output.html} \
        -j {output.json} \
        --thread {threads} 2> {log}
      touch {output.status}
    """

rule multiqc:
    input:
        # make MultiQC wait for all per-sample fastp reports to exist
        expand(OUTDIR / "{sample}" / "status" / "fastp_qc.{sample}.txt", sample=SAMPLES)
    output:
        report = OUTDIR / "multiqc" / "multiqc_report.html",
    params:
        outdir = OUTDIR / "multiqc",
        search = OUTDIR,  # Search from root to find all */qc/ subdirs
    log: OUTDIR / "log" / "multiqc.log"
    # if multiqc isn't in trim.yaml, use an env that has it (e.g., ../envs/qc.yaml)
    conda: "../envs/misc.yaml"
    shell: r"""
      mkdir -p "{params.outdir}" "$(dirname {log})"
      multiqc "{params.search}" -o "{params.outdir}" -n "multiqc_report" 2> "{log}"
    """
