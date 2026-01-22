rule denovo_index:
	message: "Indexing references"
	input:
		refset24 = REFSET,
	output:
		status = OUTDIR /"status" / "denovo.index.status"
	params:
		prefix_ref24 =  REFSET.stem,
		loc_ref24 = REFSET.parent,
	threads: config["threads"]["bowtieindex"]
	conda: "../envs/bowtie.yaml"
	shell:"""
	bowtie2-build --threads {threads} {input.refset24} {params.loc_ref24}/{params.prefix_ref24} &> /dev/null

	touch {output.status}
	"""

rule denovo_bowtie_ref24:
	message: "Align against reference set"
	input:
		r1 = rules.extract_chlamydiales_reads.output.r1,
		r2 = rules.extract_chlamydiales_reads.output.r2,
		status = rules.denovo_index.output.status,
	output:
		bam = OUTDIR / "{sample}" / "denovo" / "bowtie_ref24" / "{sample}.ref24.bam",
		coverage = OUTDIR / "{sample}" / "denovo" / "bowtie_ref24" / "{sample}.coverage.ref24.tsv",
		flagstat = OUTDIR / "{sample}" / "denovo" / "bowtie_ref24" / "{sample}.flagstat.txt",
		status = OUTDIR / "status" / "denovo.bowtie2.ref24.{sample}.txt",
	params:
		prefix = rules.denovo_index.params.prefix_ref24,
		loc = rules.denovo_index.params.loc_ref24,
		sample_name = lambda w: w.sample,
		tmp = lambda w: OUTDIR / f"{w.sample}" / "denovo" / "bowtie_ref24" / f"{w.sample}.tmp.cov.tsv",
	threads: config["threads"]["bowtie"]
	log: OUTDIR / "{sample}" / "log" / "bowtie2.ref24.{sample}.log"
	conda: "../envs/bowtie.yaml"
	shell:"""
	bowtie2 -x \
	{params.loc}/{params.prefix} \
	-1 {input.r1} \
	-2 {input.r2} \
	-p {threads} 2>> {log} | \
	samtools view -bSh - | \
	samtools sort -@{threads} \
	-o {output.bam} 2>> {log}

	samtools index {output.bam}
	samtools coverage {output.bam} > {params.tmp}

	csvtk mutate2 -t -C $ {params.tmp} \
	 -n sample_name --at 1 \
	 -e "'{params.sample_name}'" > {output.coverage}

	# Generate alignment statistics
	samtools flagstat {output.bam} > {output.flagstat}

	touch {output.status}
	"""

rule shovill:
	message: "Assembling Sample {wildcards.sample}"
	input:
		r1 = rules.extract_chlamydiales_reads.output.r1,
		r2 = rules.extract_chlamydiales_reads.output.r2,
	output:
		status = OUTDIR / "status" / "denovo.shovill.{sample}.txt",
		outdir = directory(OUTDIR / "{sample}" / "denovo" / "shovill"),
		contig = OUTDIR / "{sample}" / "denovo" / "shovill" / "contigs.fa",
	params:
		gsize = config['shovill']['gsize'],
		depth = config['shovill']['downsample'],
	conda: "../envs/shovill.yaml"
	log: OUTDIR / "{sample}" / "log" / "denovo.shovill.{sample}.log"
	benchmark: OUTDIR / "{sample}" / "benchmark" / "shovill.{sample}.txt"
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

rule scaffold:
        input:
                contigs = rules.shovill.output.contig
        output:
                outdir = directory(OUTDIR / "{sample}" / "denovo" / "scaffold"),
                scaffold = OUTDIR / "{sample}" / "denovo" / "scaffold" / "ragtag.scaffold.fasta",
                status = OUTDIR / "status" / "denovo.scaffold.{sample}.txt",
        params:
                ref = SCAFFOLDREF,
                min_len = config['ragtag']['min_len']
        threads: config['threads']['ragtag']
        conda: "../envs/scaffold.yaml"
        log: OUTDIR / "{sample}" / "log" / "denovo.scaffold.{sample}.log"
        benchmark: OUTDIR / "{sample}" / "benchmark" / "scaffold.{sample}.txt"
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

rule scaffold_summary:
    input:
        expand(OUTDIR / "{sample}" / "denovo" / "scaffold" / "ragtag.scaffold.fasta", sample=SAMPLES)
    output:
        report = OUTDIR / "reports" / "scaffold_summary.txt"
    shell: r"""
      mkdir -p "$(dirname {output.report})"
      
      echo "RagTag Scaffolding Summary Report" > {output.report}
      echo "Generated: $(date)" >> {output.report}
      echo "================================" >> {output.report}
      echo "" >> {output.report}
      
      # Count successes and failures
      success=$(find {OUTDIR}/*/denovo/scaffold -name "SCAFFOLD_SUCCESS.txt" 2>/dev/null | wc -l)
      failed=$(find {OUTDIR}/*/denovo/scaffold -name "SCAFFOLD_FAILED.txt" 2>/dev/null | wc -l)
      
      echo "Total samples processed: $(($success + $failed))" >> {output.report}
      echo "Successfully scaffolded: $success" >> {output.report}
      echo "Failed scaffolding: $failed" >> {output.report}
      echo "" >> {output.report}
      
      if [ $failed -gt 0 ]; then
        echo "Samples that FAILED scaffolding (using original contigs):" >> {output.report}
        find {OUTDIR}/*/denovo/scaffold -name "SCAFFOLD_FAILED.txt" 2>/dev/null | \
          sed 's|{OUTDIR}/||' | sed 's|/denovo/scaffold/SCAFFOLD_FAILED.txt||' | \
          sed 's/^/  - /' >> {output.report}
        echo "" >> {output.report}
      fi
      
      if [ $success -gt 0 ]; then
        echo "Samples that SUCCEEDED scaffolding:" >> {output.report}
        find {OUTDIR}/*/denovo/scaffold -name "SCAFFOLD_SUCCESS.txt" 2>/dev/null | \
          sed 's|{OUTDIR}/||' | sed 's|/denovo/scaffold/SCAFFOLD_SUCCESS.txt||' | \
          sed 's/^/  - /' >> {output.report}
      fi
      
      echo "" >> {output.report}
      echo "Check individual log files for details:" >> {output.report}
      echo "  {OUTDIR}/<sample>/log/denovo.scaffold.<sample>.log" >> {output.report}
    """

rule gapfiller:
	input:
		scaffold = rules.scaffold.output.scaffold,
		r1 = rules.extract_chlamydiales_reads.output.r1,
		r2 = rules.extract_chlamydiales_reads.output.r2,
	output:
		filled = OUTDIR / "{sample}" / "denovo" / "gap2seq" / "filled.fasta",
		status = OUTDIR / "status" / "denovo.gap2seq.{sample}.txt",
	shadow: "minimal"
	threads: config['threads']['gap2seq']
	conda: "../envs/scaffold.yaml"
	log: OUTDIR / "{sample}" / "log" / "denovo.gap2seq.{sample}.log"
	benchmark: OUTDIR / "{sample}" / "benchmark" / "denovo.gap2seq.{sample}.txt"
	shell:"""
	set +e
	Gap2Seq \
	--scaffolds {input.scaffold} \
	--filled {output.filled} \
	--reads {input.r1},{input.r2} \
	--threads {threads}	> {log} 2>&1
	exitcode=$?

	# gap2seq will fail if filling fails.
	# capture error, copy scaffolded as filled.
	if [ $exitcode -ne 0 ]
	then
		cp {input.scaffold} {output.filled}
	fi

	touch {output.status}
	"""

rule rename:
	input:
		filled = rules.gapfiller.output.filled
	output:
		renamed = OUTDIR / "{sample}" / "denovo" / "denovo_{sample}.final.fasta",
	params:
		sample_name = lambda w: w.sample
	conda: "../envs/misc.yaml"
	shell:"""
	seqkit replace -p "(.*)" -r '{params.sample_name}_contig{{nr}}' {input.filled} > {output.renamed} 
	"""

rule quast:
    message: "genome assembly statistics"
    input:
        contig  = rules.rename.output.renamed
    output:
        quasted = directory(OUTDIR / "{sample}" / "denovo" / "assembly_statistics"),
        report_tsv = OUTDIR / "{sample}" / "denovo" / "assembly_statistics" / "report.tsv",  # ← ADDED
        status = OUTDIR / "status" / "denovo.assembly_statistics.{sample}.txt",
    conda: "../envs/misc.yaml"
    log: OUTDIR / "{sample}" / "log" / "denovo.assembly_statistics.{sample}.log"
    threads: config["threads"]["quast"]
    shell:"""
    quast \
    {input.contig} \
    -t {threads} \
    -o {output.quasted} > {log} 2>&1

    touch {output.status}
    """

rule denovo_collate_coverage:
	input:
		coverages = expand(OUTDIR / "{sample}" / "denovo" / "bowtie_ref24" / "{sample}.coverage.ref24.tsv", sample = SAMPLES),
	output:
		coverages = OUTDIR / "denovo.coverage.tsv",
		status = OUTDIR / "status" / "denovo.collate.coverage.txt",
	threads: 1
	conda: "../envs/misc.yaml"
	shell:"""
	csvtk concat  {input.coverages} > {output.coverages}

	touch {output.status}
	"""

rule denovo_collate_flagstat:
	input:
		flagstats = expand(OUTDIR / "{sample}" / "denovo" / "bowtie_ref24" / "{sample}.flagstat.txt", sample = SAMPLES),
	output:
		summary = OUTDIR / "reports" / "denovo.alignment_summary.txt",
		status = OUTDIR / "status" / "denovo.collate.flagstat.txt",
	threads: 1
	shell: r"""
	mkdir -p "$(dirname {output.summary})"

	echo "De Novo Alignment Summary (vs Reference Set)" > {output.summary}
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

