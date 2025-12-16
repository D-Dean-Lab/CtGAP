#include: "/home/shola/CtGAP/workflow/rules/5-align.smk"
#include: "/home/shola/CtGAP/workflow/rules/4-assemble.smk"

rule dna_diff:
	message: "dnadiff comparison between denovo and reference-guided assemblies"
	input: 
		denovo = OUTDIR / "{sample}" / "denovo" / "denovo_{sample}.final.fasta",
		ref_denovo = OUTDIR / "{sample}" / "ref-denovo" / "ref-denovo_{sample}.final.fasta",
	output: 
		dna_diff = OUTDIR / "{sample}" / "dnadiff_report" / "{sample}.report",
		status = OUTDIR / "status" / "dnadiff.{sample}.txt",
	conda: "../envs/misc.yaml"
	log: OUTDIR / "{sample}" / "log" / "dnadiff.{sample}.log"
	threads: 5
	params:
		prefix = directory(OUTDIR / "{sample}" / "dnadiff_report" / "{sample}"),
	shell:"""
	mkdir -p "$(dirname {params.prefix})"
	
	perl $CONDA_PREFIX/bin/dnadiff \
	{input.denovo} \
	{input.ref_denovo} \
	-p {params.prefix} > {log} 2>&1

	touch {output.status}
	"""
