include: "4-assemble.smk"

rule quast:
	message: "genome quality check"
	input:
		fasta = rules.rename.output.renamed,
	output:
		outdir = directory(OUTDIR/"{sample}"/"quast"/"quast_output.html")
		#status = OUTDIR / "status" / "denovo.quasted.{sample}.txt",
	conda: "../envs/misc.yaml"
	log: OUTDIR / "{sample}" / "log" / "quast.{sample}.log"
	threads: config["threads"]["quast"]
	shell: """
	quast \
	{input.fasta} \
	-t {threads} \
	-o {output.outdir}

	touch {output.status}
	"""
