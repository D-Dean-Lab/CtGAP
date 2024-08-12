include: "/home/oolago2/assembly_pipeline/fullTest_Dec23/CtGAP_test/ctgap/workflow/rules/5-align.smk"
include: "/home/oolago2/assembly_pipeline/fullTest_Dec23/CtGAP_test/ctgap/workflow/rules/4-assemble.smk"

rule dna_diff:
	message: "dnadiff comparison between denovo and reference denovo assemblies"
	input: 
		denovo = rules.rename.output.renamed,
		ref_denovo = rules.ref_rename.output.renamed,
	output: 
		dna_diff = OUTDIR / "{sample}" / "dnadiff_report" / "{sample}.report",
		status = OUTDIR / "status" / "dnadiff.{sample}.txt",
	conda: "../envs/misc.yaml"
	log: OUTDIR / "{sample}" / "log" / "dnadiff.{sample}.log"
	threads: 5
	params:
		prefix = directory(OUTDIR / "{sample}" / "dnadiff_report" / "{sample}"),
	shell:"""
	perl $CONDA_PREFIX/bin/dnadiff \
	{input.denovo} \
	{input.ref_denovo} \
	-p {params.prefix} > {log} 2>&1


	touch {output.status}
	"""