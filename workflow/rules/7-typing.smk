# workflow/rules/7-typing.smk
# Typing and downstream analyses - works with any assembly mode

def get_best_assembly(wildcards):
    """Return the best assembly path based on mode"""
    return OUTDIR / wildcards.sample / "best" / "assembly.fasta"

rule blast_ompa:
    input:
        contig = get_best_assembly,  # ← Uses best assembly
    output:
        primary_tab = OUTDIR / "{sample}" / "best" / "blast" / "blast.ompa.tab",
        secondary_tab = OUTDIR / "{sample}" / "best" / "blast" / "blast.secondary.tab",
        status = OUTDIR / "status" / "best.blastn.{sample}.txt"
    log: OUTDIR / "{sample}" / "log" / "best.blast.ompa.{sample}.log"
    params:
        outfmt = 6,
        primary_db = OMPABLASTDB,
        secondary_db = config.get("secondaryBlastDB", "resources/references/ct/secondaryBlastDB/secondary"),
        targets = 1,
        target_hits = "A_DQ064279.1 B_M33636.1 Ba_DQ064282.1 C_AF352789.1"
    threads: config["threads"]["blast"]
    conda: "../envs/misc.yaml"
    shell: r"""
    mkdir -p "$(dirname {output.primary_tab})"
    
    blastn \
      -query {input.contig} \
      -db {params.primary_db} \
      -max_target_seqs {params.targets} \
      -outfmt {params.outfmt} \
      -out {output.primary_tab} \
      > {log} 2>&1

    # Get top hit
    top_hit=$(head -n1 {output.primary_tab} | cut -f2)
    
    # Check if top hit matches targets and run secondary BLAST
    if echo "{params.target_hits}" | grep -qw "$top_hit"; then
      echo "Top hit $top_hit matched - running secondary BLAST" >> {log}
      blastn \
        -query {input.contig} \
        -db {params.secondary_db} \
        -max_target_seqs 5 \
        -outfmt {params.outfmt} \
        -out {output.secondary_tab} \
        2>> {log}
    else
      echo "Top hit $top_hit - no secondary BLAST needed" >> {log}
      touch {output.secondary_tab}
    fi

    touch {output.status}
    """

rule mlst:
    input:
        contig = get_best_assembly,  # ← Uses best assembly
    output:
        generic = OUTDIR / "{sample}" / "best" / "mlst" / "{sample}.genome.chlamydiales.mlst.txt",
        ct = OUTDIR / "{sample}" / "best" / "mlst" / "{sample}.genome.ctrachomatis.mlst.txt",
        status = OUTDIR / "status" / "best.mlst.{sample}.txt",
    log: OUTDIR / "{sample}" / "log" / "best.mlst.{sample}.log"
    benchmark: OUTDIR / "{sample}" / "benchmark" / "best.mlst.{sample}.txt"
    conda: "../envs/mlst.yaml"
    params:
        dbgeneric = MLSTDB / "chlamydiales",
        dbct = MLSTDB / "c.trachomatis",
    threads: config["threads"]["mlst"]
    shell:"""
    echo -e "chlamydiales\n" >> {log}
    claMLST search \
    {params.dbgeneric} \
    {input} > {output.generic} 2>> {log}

    echo -e "\nctrachomatis\n" >> {log}
    claMLST search \
    {params.dbct} \
    {input} > {output.ct} 2>> {log}

    touch {output.status}
    """

rule collate_blast:
    input:
        status_blast = expand(OUTDIR / "status" / "best.blastn.{sample}.txt", sample = SAMPLES),
    output:
        tsv = OUTDIR / "best.ompA_genovar.blast.tsv",
        status = OUTDIR / "status" / "best.ompA_genovar.collate.blast.txt",
    params:
        outdir = OUTDIR,
        pattern = "**/best/blast/*.tab",
    threads: 1
    shell:"""
    echo -e "query\tsubject\tpident\tlength\tmismatch\tgapopen\tquery_start\tquery_end\tsubject_start\tsubject_end\tevalue\tbitscore" > {output.tsv}
    grep "" {params.outdir}/{params.pattern} >> {output.tsv}

    touch {output.status}
    """

rule collate_secondary_blast:
    input:
        status_blast = expand(OUTDIR / "status" / "best.blastn.{sample}.txt", sample = SAMPLES),
    output:
        tsv = OUTDIR / "best.secondary.blast.tsv",
        status = OUTDIR / "status" / "best.secondary.collate.blast.txt",
    params:
        outdir = OUTDIR,
        pattern = "**/best/blast/blast.secondary.tab",
    threads: 1
    shell: r"""
    echo -e "query\tsubject\tpident\tlength\tmismatch\tgapopen\tquery_start\tquery_end\tsubject_start\tsubject_end\tevalue\tbitscore" > {output.tsv}
    
    # Only add non-empty files
    find {params.outdir} -path "{params.pattern}" -type f ! -empty -exec cat {{}} \; >> {output.tsv}

    touch {output.status}
    """

rule collate_mlst:
    input:
        generic = expand(OUTDIR / "{sample}" / "best" / "mlst" / "{sample}.genome.chlamydiales.mlst.txt", sample = SAMPLES),
        ct = expand(OUTDIR / "{sample}" / "best" / "mlst" / "{sample}.genome.ctrachomatis.mlst.txt", sample = SAMPLES),
    output:
        generic = OUTDIR / "best.mlst.generic.results.tsv",
        ct = OUTDIR / "best.mlst.ct.results.tsv",
        status = OUTDIR / "status" / "best.mlst.collate.txt"
    conda: "../envs/misc.yaml"
    threads: 1
    shell:"""
    csvtk concat {input.generic} -o {output.generic}
    csvtk concat {input.ct} -o {output.ct}

    touch {output.status}
    """

rule ska_input_prep:
    input:
        assemblies = expand(OUTDIR / "{sample}" / "best" / "assembly.fasta", sample = SAMPLES),  # ← Changed
        individual = INDIVIDUAL,
    output:
        glist = OUTDIR / "tree" / "input.list"
    params:
        samples = SAMPLES,
    threads: 10
    run:
        for f,s in zip(input.assemblies + list(INDIVIDUAL_DICT.values()), list(params.samples) + list(INDIVIDUAL_DICT.keys())):
            with open(output.glist, "a") as handle:
                handle.write(f"{s}\t{f}\n")

rule ska_alignment:
	input:
		glist = rules.ska_input_prep.output.glist,
	output:
		alignment = OUTDIR / "tree" / "ska_alignment.fasta",
	params:
		reference = SCAFFOLDREF,
		organism = ORG,
	conda: "../envs/tree.yaml"
	log: OUTDIR / "tree" / "ska.log"
	threads: config["threads"]["ska"]
	shell:"""
	workflow/scripts/generate_ska_alignment.py \
	--input {input.glist} \
	--out {output.alignment} \
	--reference {params.reference} \
	--threads {threads} > {log} 2>&1
	"""

rule tree:
	input:
		alignment = rules.ska_alignment.output.alignment
	output:
		tree = OUTDIR / f"{ORG}.tree"
	params:
		replicates = 1000,
		prefix = ORG
	conda: "../envs/tree.yaml"
	threads: config["threads"]["tree"]
	log: OUTDIR / "tree" / "iqtree.log"
	shell:"""
	augur tree \
	--alignment {input.alignment} \
	--output {output.tree} \
	--override-default-args \
	--tree-builder-args="-alrt {params.replicates} -B {params.replicates}" \
	--nthreads {threads} > {log} 2>&1
	"""
