# workflow/rules/5-align.smk
# Reference-guided assembly pipeline 

rule index:
        message: "Indexing references"
        input:
                reference = REF,
                refset24 = REFSET,
        output:
                status = OUTDIR /"status" / "ref-denovo.index.status"
        params:
                prefix_ref = REF.stem,
                prefix_ref24 =  REFSET.stem,
                loc_ref = REF.parent,
                loc_ref24 = REFSET.parent,
        threads: config["threads"]["bowtieindex"]
        conda: "../envs/bowtie.yaml"
        shell:"""
        bowtie2-build --threads {threads} {input.reference} {params.loc_ref}/{params.prefix_ref} &> /dev/null
        bowtie2-build --threads {threads} {input.refset24} {params.loc_ref24}/{params.prefix_ref24} &> /dev/null

        touch {output.status}
        """

rule bowtie:
        message: "Aligning to specified reference"
        input:
                r1 = rules.extract_chlamydiales_reads.output.r1,
                r2 = rules.extract_chlamydiales_reads.output.r2,
                status = rules.index.output.status,
        output:
                bam = OUTDIR / "{sample}" / "ref-denovo" / "bowtie2" / "{sample}.bam",
                flagstat = OUTDIR / "{sample}" / "ref-denovo" / "bowtie2" / "{sample}.flagstat.txt",
                status = OUTDIR / "status" / "ref-denovo.bowtie2.{sample}.txt",
        params:
                prefix = rules.index.params.prefix_ref,
                loc = rules.index.params.loc_ref,
        threads: config["threads"]["bowtie"]
        log: OUTDIR / "{sample}" / "log" / "ref-denovo.bowtie2.{sample}.log"
        conda: "../envs/bowtie.yaml"
        shell:"""
        bowtie2 -x \
        {params.loc}/{params.prefix} \
        -1 {input.r1} \
        -2 {input.r2} \
        --threads {threads} 2>> {log} | \
        samtools view -bSh - | \
        samtools sort -n -@{threads} \
        -o {output.bam} 2>> {log}

        # Generate alignment statistics
        samtools flagstat {output.bam} > {output.flagstat}

        touch {output.status}
        """

rule bowtie_ref24:
        message: "Align against reference set"
        input:
                r1 = rules.extract_chlamydiales_reads.output.r1,
                r2 = rules.extract_chlamydiales_reads.output.r2,
                status = rules.index.output.status,
        output:
                bam = OUTDIR / "{sample}" / "ref-denovo" / "bowtie_ref24" / "{sample}.ref24.bam",
                coverage = OUTDIR / "{sample}" / "ref-denovo" / "bowtie_ref24" / "{sample}.coverage.ref24.tsv",
                status = OUTDIR / "status" / "ref-denovo.bowtie2.ref24.{sample}.txt",
        params:
                prefix = rules.index.params.prefix_ref24,
                loc = rules.index.params.loc_ref24,
                sample_name = lambda w: w.sample,
                tmp = lambda w: OUTDIR / f"{w.sample}" / "ref-denovo" / "bowtie_ref24" / f"{w.sample}.tmp.cov.tsv",
        threads: config["threads"]["bowtie"]
        log: OUTDIR / "{sample}" / "log" / "ref-denovo.bowtie2.ref24.{sample}.log"
        conda: "../envs/bowtie.yaml"
        shell:"""
        bowtie2 -x \
        {params.loc}/{params.prefix} \
        -1 {input.r1} \
        -2 {input.r2} \
        --threads {threads} 2>> {log} | \
        samtools view -bSh - | \
        samtools sort -@{threads} \
        -o {output.bam} 2>> {log}

        samtools index {output.bam}
        samtools coverage {output.bam} > {params.tmp}

        csvtk mutate2 -t -C $ {params.tmp} -n "sample_name" --at 1 -e " '{params.sample_name}' " > {output.coverage}

        touch {output.status}
        """

rule bamtofastq:
        message: "Extracting aligned reads from reference-guided alignment"
        input:
                bam = rules.bowtie.output.bam,
        output:
                r1 = OUTDIR / "{sample}" / "ref-denovo" / "bamtofastq" / "{sample}_refg_r1.fastq.gz",
                r2 = OUTDIR / "{sample}" / "ref-denovo" / "bamtofastq" / "{sample}_refg_r2.fastq.gz",
                status = OUTDIR / "status" / "ref-denovo.bamtofastq.{sample}.txt",
        params:
                r1_tmp = lambda w: str(OUTDIR / w.sample / "ref-denovo" / "bamtofastq" / f"{w.sample}_refg_r1.fastq"),
                r2_tmp = lambda w: str(OUTDIR / w.sample / "ref-denovo" / "bamtofastq" / f"{w.sample}_refg_r2.fastq"),
        threads: config["threads"]["bedtools"]
        log: OUTDIR / "{sample}" / "log" / "ref-denovo.bamtofastq.{sample}.log"
        conda: "../envs/misc.yaml"
        shell:"""
        bedtools bamtofastq \
        -i {input.bam} \
        -fq {params.r1_tmp} \
        -fq2 {params.r2_tmp} 2> {log}

        gzip -f {params.r1_tmp}
        gzip -f {params.r2_tmp}

        touch {output.status}
        """

rule ref_shovill:
        message: "Assembling reference-guided reads for sample {wildcards.sample}"
        input:
                r1 = rules.bamtofastq.output.r1,
                r2 = rules.bamtofastq.output.r2,
        output:
                status = OUTDIR / "status" / "ref-denovo.shovill.{sample}.txt",
                outdir = directory(OUTDIR / "{sample}" / "ref-denovo" / "shovill"),
                contig = OUTDIR / "{sample}" / "ref-denovo" / "shovill" / "contigs.fa",
        params:
                gsize = config['shovill']['gsize'],
                depth = config['shovill']['downsample'],
        conda: "../envs/shovill.yaml"
        log: OUTDIR / "{sample}" / "log" / "ref-denovo.shovill.{sample}.log"
        benchmark: OUTDIR / "{sample}" / "benchmark" / "ref-denovo.shovill.{sample}.txt"
        threads: config["threads"]["shovill"]
        shell:"""
        shovill \
        --R1 {input.r1} \
        --R2 {input.r2} \
        --outdir {output.outdir} \
        --gsize {params.gsize} \
        --depth {params.depth} \
        --force \
        --cpus {threads} 2>> {log} 1>> {log}

        touch {output.status}
        """

rule ref_scaffold:
        message: "Scaffolding reference-guided assembly against reference"
        input:
                contigs = rules.ref_shovill.output.contig
        output:
                outdir = directory(OUTDIR / "{sample}" / "ref-denovo" / "scaffold"),
                scaffold = OUTDIR / "{sample}" / "ref-denovo" / "scaffold" / "ragtag.scaffold.fasta",
                status = OUTDIR / "status" / "ref-denovo.scaffold.{sample}.txt",
        params:
                ref = SCAFFOLDREF,
                min_len = config['ragtag']['min_len']
        threads: config['threads']['ragtag']
        conda: "../envs/scaffold.yaml"
        log: OUTDIR / "{sample}" / "log" / "ref-denovo.scaffold.{sample}.log"
        benchmark: OUTDIR / "{sample}" / "benchmark" / "ref-denovo.scaffold.{sample}.txt"
        shell: r"""
          ragtag.py scaffold \
            {params.ref} \
            {input.contigs} \
            -o {output.outdir} \
            -f {params.min_len} \
            -r \
            -u \
            -t {threads} \
            2> {log} || true

          # If ragtag failed, copy original contigs as scaffold
          if [ ! -f {output.scaffold} ]; then
            echo "WARNING: RagTag scaffolding failed for {wildcards.sample}. Using original contigs." >> {log}
            echo "SCAFFOLD_FAILED" > {output.outdir}/SCAFFOLD_FAILED.txt
            cp {input.contigs} {output.scaffold}
          else
            echo "SCAFFOLD_SUCCESS" > {output.outdir}/SCAFFOLD_SUCCESS.txt
          fi

          touch {output.status}
        """

rule ref_scaffold_summary:
    input:
        expand(OUTDIR / "{sample}" / "ref-denovo" / "scaffold" / "ragtag.scaffold.fasta", sample=SAMPLES)
    output:
        report = OUTDIR / "reports" / "ref_scaffold_summary.txt"
    shell: r"""
      mkdir -p "$(dirname {output.report})"

      echo "RagTag Scaffolding Summary Report (Reference-Guided Mode)" > {output.report}
      echo "Generated: $(date)" >> {output.report}
      echo "================================" >> {output.report}
      echo "" >> {output.report}

      # Count successes and failures
      success=$(find {OUTDIR}/*/ref-denovo/scaffold -name "SCAFFOLD_SUCCESS.txt" 2>/dev/null | wc -l)
      failed=$(find {OUTDIR}/*/ref-denovo/scaffold -name "SCAFFOLD_FAILED.txt" 2>/dev/null | wc -l)

      echo "Total samples processed: $(($success + $failed))" >> {output.report}
      echo "Successfully scaffolded: $success" >> {output.report}
      echo "Failed scaffolding: $failed" >> {output.report}
      echo "" >> {output.report}

      if [ $failed -gt 0 ]; then
        echo "Samples that FAILED scaffolding (using original contigs):" >> {output.report}
        find {OUTDIR}/*/ref-denovo/scaffold -name "SCAFFOLD_FAILED.txt" 2>/dev/null | \
          sed 's|{OUTDIR}/||' | sed 's|/ref-denovo/scaffold/SCAFFOLD_FAILED.txt||' | \
          sed 's/^/  - /' >> {output.report}
        echo "" >> {output.report}
      fi

      if [ $success -gt 0 ]; then
        echo "Samples that SUCCEEDED scaffolding:" >> {output.report}
        find {OUTDIR}/*/ref-denovo/scaffold -name "SCAFFOLD_SUCCESS.txt" 2>/dev/null | \
          sed 's|{OUTDIR}/||' | sed 's|/ref-denovo/scaffold/SCAFFOLD_SUCCESS.txt||' | \
          sed 's/^/  - /' >> {output.report}
      fi

      echo "" >> {output.report}
      echo "Check individual log files for details:" >> {output.report}
      echo "  {OUTDIR}/<sample>/log/ref-denovo.scaffold.<sample>.log" >> {output.report}
    """

rule ref_gapfiller:
        message: "Gap filling reference-guided assembly"
        input:
                scaffold = rules.ref_scaffold.output.scaffold,
                r1 = rules.extract_chlamydiales_reads.output.r1,
                r2 = rules.extract_chlamydiales_reads.output.r2,
        output:
                filled = OUTDIR / "{sample}" / "ref-denovo" / "gap2seq" / "{sample}.fasta",
                status = OUTDIR / "status" / "ref-denovo.gap2seq.{sample}.txt",
        shadow: "minimal"
        threads: config['threads']['gap2seq']
        conda: "../envs/scaffold.yaml"
        log: OUTDIR / "{sample}" / "log" / "ref-denovo.gap2seq.{sample}.log"
        benchmark: OUTDIR / "{sample}" / "benchmark" / "ref-denovo.gap2seq.{sample}.txt"
        shell:"""
        set +e
        Gap2Seq \
        --scaffolds {input.scaffold} \
        --filled {output.filled} \
        --reads {input.r1},{input.r2} \
        --threads {threads} > {log} 2>&1
        exitcode=$?

        # gap2seq will fail if filling fails.
        # capture error, copy scaffolded as filled.
        if [ $exitcode -ne 0 ]
        then
                echo "Gap2Seq failed, copying scaffold as filled sequence" >> {log}
                cp {input.scaffold} {output.filled}
        fi

        touch {output.status}
        """

rule ref_rename:
        message: "Renaming contigs for reference-guided assembly"
        input:
                filled = rules.ref_gapfiller.output.filled
        output:
                renamed = OUTDIR / "{sample}" / "ref-denovo" / "ref-denovo_{sample}.final.fasta",
        params:
                sample_name = lambda w: w.sample
        conda: "../envs/misc.yaml"
        shell:"""
        seqkit replace -p "(.*)" -r '{params.sample_name}_contig{{nr}}' {input.filled} > {output.renamed}
        """

rule ref_quast:
    message: "Genome assembly statistics for reference-guided assembly (reference-free)"
    input:
        contig  = rules.ref_rename.output.renamed,
    output:
        quasted = directory(OUTDIR / "{sample}" / "ref-denovo" / "quast"),
        report_tsv = OUTDIR / "{sample}" / "ref-denovo" / "quast" / "report.tsv",  # ← ADDED
        status = OUTDIR / "status" / "ref-denovo.assembly_statistics.{sample}.txt",
    conda: "../envs/misc.yaml"
    log: OUTDIR / "{sample}" / "log" / "ref-denovo.assembly_statistics.{sample}.log"
    threads: config["threads"]["quast"]
    shell:"""
    quast \
    {input.contig} \
    -t {threads} \
    -o {output.quasted} > {log} 2>&1

    touch {output.status}
    """

rule ref_quast_scaffolded:
        message: "Genome assembly statistics vs scaffold reference"
        input:
                contig  = rules.ref_rename.output.renamed,
        output:
                quasted = directory(OUTDIR / "{sample}" / "ref-denovo" / "assembly_statistics_vs_scaffold"),
                status = OUTDIR / "status" / "ref-denovo.assembly_statistics_vs_scaffold.{sample}.txt",
        params:
                reference = SCAFFOLDREF,
        conda: "../envs/misc.yaml"
        log: OUTDIR / "{sample}" / "log" / "ref-denovo.assembly_statistics_vs_scaffold.{sample}.log"
        threads: config["threads"]["quast"]
        shell:"""
        quast \
        {input.contig} \
        -R {params.reference} \
        -t {threads} \
        -o {output.quasted} > {log} 2>&1

        touch {output.status}
        """

rule ref_collate_coverage:
        input:
                coverages = expand(OUTDIR / "{sample}" / "ref-denovo" / "bowtie_ref24" / "{sample}.coverage.ref24.tsv", sample = SAMPLES),
        output:
                coverages = OUTDIR / "ref-denovo.coverage.tsv",
                status = OUTDIR / "status" / "ref-denovo.collate.coverage.txt",
        threads: 1
        conda: "../envs/misc.yaml"
        shell:"""
        csvtk concat -C $ {input.coverages} > {output.coverages}

        touch {output.status}
        """

rule ref_collate_flagstat:
        input:
                flagstats = expand(OUTDIR / "{sample}" / "ref-denovo" / "bowtie2" / "{sample}.flagstat.txt", sample = SAMPLES),
        output:
                summary = OUTDIR / "reports" / "ref-denovo.alignment_summary.txt",
                status = OUTDIR / "status" / "ref-denovo.collate.flagstat.txt",
        threads: 1
        shell: r"""
        mkdir -p "$(dirname {output.summary})"
        
        echo "Reference-Guided Alignment Summary" > {output.summary}
        echo "Generated: $(date)" >> {output.summary}
        echo "================================" >> {output.summary}
        echo "" >> {output.summary}
        
        for flagstat in {input.flagstats}; do
            sample=$(basename $(dirname $(dirname $flagstat)))
            echo "Sample: $sample" >> {output.summary}
            cat $flagstat >> {output.summary}
            echo "" >> {output.summary}
        done
        
        touch {output.status}
        """
