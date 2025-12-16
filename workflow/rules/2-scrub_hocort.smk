# workflow/rules/2-scrub_hocort.smk
import os

# --- 1) Build the host index *with HoCoRT* (once) ---
rule index_host_hocort:
    input:
        fasta = HOSTREF
    output:
        idx = [
            "resources/grch38.bt2/grch38.1.bt2",
            "resources/grch38.bt2/grch38.2.bt2",
            "resources/grch38.bt2/grch38.3.bt2",
            "resources/grch38.bt2/grch38.4.bt2",
            "resources/grch38.bt2/grch38.rev.1.bt2",
            "resources/grch38.bt2/grch38.rev.2.bt2",
        ]
    params:
        outprefix = "resources/grch38.bt2/grch38"
    threads: 4
    conda: "../envs/hocort.yaml"
    shell: r"""
      mkdir -p "$(dirname {params.outprefix})"
      hocort index bowtie2 --input {input.fasta} --output {params.outprefix}
    """

# Deplete host reads with HoCoRT (bowtie2 backend)
rule host_deplete_hocort:
    input:
        r1  = rules.fastp.output.r1,
        r2  = rules.fastp.output.r2,
        idx = rules.index_host_hocort.output.idx
    output:
        r1     = OUTDIR / "{sample}" / "scrub" / "{sample}_hostdepleted_R1.fastq.gz",
        r2     = OUTDIR / "{sample}" / "scrub" / "{sample}_hostdepleted_R2.fastq.gz",
        status = OUTDIR / "status" / "hocort.{sample}.txt"
    params:
        prefix        = "resources/grch38.bt2/grch38",
        opts          = config.get("hocort_opts", ""),
        bowtie2_args  = config.get("hocort_bowtie2_args", "")
    threads: config["threads"]["hocort"]
    conda: "../envs/hocort.yaml"
    log: OUTDIR / "{sample}" / "log" / "hocort.{sample}.log"
    benchmark: OUTDIR / "{sample}" / "benchmark" / "hocort.{sample}.txt"
    shell: r"""
      set -euo pipefail
      mkdir -p "$(dirname {output.r1})" "$(dirname {log})" "$(dirname {output.status})"

      # Build optional -c flag only if args are provided (avoids commas error)
      EXTRA_C=""
      if [ -n "{params.bowtie2_args}" ]; then
        EXTRA_C='-c="{params.bowtie2_args}"'
      fi

      # HoCoRT host depletion: keep UNMAPPED reads (-f true), end-to-end mode
      hocort map bowtie2 \
        -x {params.prefix} \
        -i {input.r1} {input.r2} \
        -o {output.r1} {output.r2} \
        -t {threads} \
        -p end-to-end \
        -f true \
        $EXTRA_C \
        {params.opts} \
        &> {log}

      touch {output.status}
    """


# Classify host-depleted reads with Kraken2
rule kraken2_classify:
    input:
        r1 = rules.host_deplete_hocort.output.r1,
        r2 = rules.host_deplete_hocort.output.r2
    output:
        report = OUTDIR / "{sample}" / "scrub" / "{sample}.kraken2.report",
        labels = OUTDIR / "{sample}" / "scrub" / "{sample}.kraken2.labels",
        k2 = OUTDIR / "{sample}" / "scrub" / "{sample}.kraken2"
    params:
        db = KRAKENDB
    threads: config["threads"]["kraken"]
    conda: "../envs/kraken.yaml"
    log: OUTDIR / "{sample}" / "log" / "kraken2.{sample}.log"
    benchmark: OUTDIR / "{sample}" / "benchmark" / "kraken2.{sample}.txt"
    shell: r"""
      mkdir -p "$(dirname {output.report})" "$(dirname {log})"

      kraken2 \
        --db {params.db} \
        --threads {threads} \
        --report {output.report} \
        --paired {input.r1} {input.r2} \
        --output {output.k2} \
        --use-names  
        &> {log}

      # compact labels for KrakenTools (handle empty output)
      if [ -s {output.k2} ]; then
        cut -f2,3 {output.k2} > {output.labels}
      else
        touch {output.labels}
      fi
    """

# Extract Chlamydiales (+ descendants) with KrakenTools
rule extract_chlamydiales_reads:
    input:
        r1   = rules.host_deplete_hocort.output.r1,
        r2   = rules.host_deplete_hocort.output.r2,
        k2   = rules.kraken2_classify.output.k2,
        rep  = rules.kraken2_classify.output.report  # <- pass the .report explicitly
    output:
        r1     = OUTDIR / "{sample}" / "scrub" / "{sample}_trim_scrub_R1.fastq.gz",
        r2     = OUTDIR / "{sample}" / "scrub" / "{sample}_trim_scrub_R2.fastq.gz",
        status = OUTDIR / "status" / "scrub.{sample}.txt"
    params:
        taxid = config.get("chlamydiales_taxid", 51291)
    threads: 2
    conda: "../envs/kraken.yaml"
    log: OUTDIR / "{sample}" / "log" / "extract_ct.{sample}.log"
    benchmark: OUTDIR / "{sample}" / "benchmark" / "extract_ct.{sample}.txt"
    shell: r"""
      set -euo pipefail
      mkdir -p "$(dirname {output.r1})" "$(dirname {output.r2})" "$(dirname {log})" "$(dirname {output.status})"

      tmp1="$(mktemp)"; tmp2="$(mktemp)"
      extract_kraken_reads.py \
        -k {input.k2} \
        --report {input.rep} \
        -s {input.r1} -s2 {input.r2} \
        -o "$tmp1" -o2 "$tmp2" \
        -t {params.taxid} \
        --include-children \
        --include-parents \
        --fastq-output \
        &> {log}

      gzip -c "$tmp1" > {output.r1}
      gzip -c "$tmp2" > {output.r2}
      rm -f "$tmp1" "$tmp2"

      touch {output.status}
    """
