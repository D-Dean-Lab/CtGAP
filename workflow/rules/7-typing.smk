# workflow/rules/7-typing.smk
# Typing and downstream analyses - works with any assembly mode

def get_best_assembly(wildcards):
    """Return the best assembly path based on mode"""
    return OUTDIR / wildcards.sample / "best" / "assembly.fasta"

# -----------------------------------------------------------------------------
# BLAST DATABASE PREPARATION
# Ensures databases are built with the same BLAST version used for searching
# -----------------------------------------------------------------------------

rule prepare_ompa_blastdb:
    """Build/rebuild ompA BLAST database to ensure version compatibility."""
    input:
        fasta = Path(config["ompaBlastDB"]).parent / (Path(config["ompaBlastDB"]).name + ".fasta")
    output:
        status = OUTDIR / "status" / "blastdb.ompa.txt"
    params:
        db_prefix = config["ompaBlastDB"]
    log: OUTDIR / "log" / "prepare_ompa_blastdb.log"
    conda: "../envs/misc.yaml"
    shell: r"""
    # Always rebuild to ensure version compatibility with this conda env's BLAST
    echo "Building ompA BLAST database with $(blastn -version | head -1)..." > {log}
    makeblastdb -in {input.fasta} -dbtype nucl -out {params.db_prefix} >> {log} 2>&1
    echo "Database built successfully" >> {log}
    touch {output.status}
    """

rule prepare_secondary_blastdb:
    """Build/rebuild secondary BLAST database to ensure version compatibility."""
    input:
        fasta = Path(config.get("secondaryBlastDB", "resources/references/ct/secondaryBlastDB/secondary")).parent / (Path(config.get("secondaryBlastDB", "resources/references/ct/secondaryBlastDB/secondary")).name + ".fasta")
    output:
        status = OUTDIR / "status" / "blastdb.secondary.txt"
    params:
        db_prefix = config.get("secondaryBlastDB", "resources/references/ct/secondaryBlastDB/secondary")
    log: OUTDIR / "log" / "prepare_secondary_blastdb.log"
    conda: "../envs/misc.yaml"
    shell: r"""
    # Always rebuild to ensure version compatibility with this conda env's BLAST
    echo "Building secondary BLAST database with $(blastn -version | head -1)..." > {log}
    makeblastdb -in {input.fasta} -dbtype nucl -out {params.db_prefix} >> {log} 2>&1
    echo "Database built successfully" >> {log}
    touch {output.status}
    """

# -----------------------------------------------------------------------------
# BLAST TYPING
# -----------------------------------------------------------------------------

rule blast_ompa:
    input:
        contig = get_best_assembly,  # ← Uses best assembly
        # Ensure BLAST databases are built with compatible version
        ompa_db_status = OUTDIR / "status" / "blastdb.ompa.txt",
        secondary_db_status = OUTDIR / "status" / "blastdb.secondary.txt",
        # Wait for assembly filtering to complete (if enabled)
        filter_status = OUTDIR / "status" / "filter_applied.{sample}.txt",
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

    # Primary BLAST against ompA database
    blastn \
      -query {input.contig} \
      -db {params.primary_db} \
      -num_threads {threads} \
      -use_index false \
      -max_target_seqs {params.targets} \
      -outfmt {params.outfmt} \
      -out {output.primary_tab} \
      > {log} 2>&1

    # Check for empty results
    if [ ! -s {output.primary_tab} ]; then
      echo "WARNING: No BLAST hits found for {wildcards.sample}" >> {log}
      touch {output.secondary_tab}
      touch {output.status}
      exit 0
    fi

    # Get top hit by bitscore (column 12)
    top_hit=$(sort -k12,12nr {output.primary_tab} | head -n1 | cut -f2)

    # Check if top hit matches A/B/Ba/C genotypes requiring secondary BLAST
    if echo "{params.target_hits}" | grep -qw "$top_hit"; then
      echo "Top hit $top_hit matched - running secondary BLAST" >> {log}
      blastn \
        -query {input.contig} \
        -db {params.secondary_db} \
        -num_threads {threads} \
        -use_index false \
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
        # Wait for assembly filtering to complete (if enabled)
        filter_status = OUTDIR / "status" / "filter_applied.{sample}.txt",
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
    {input.contig} > {output.generic} 2>> {log}

    echo -e "\nctrachomatis\n" >> {log}
    claMLST search \
    {params.dbct} \
    {input.contig} > {output.ct} 2>> {log}

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
    threads: 1
    shell: r"""
    echo -e "query\tsubject\tpident\tlength\tmismatch\tgapopen\tquery_start\tquery_end\tsubject_start\tsubject_end\tevalue\tbitscore" > {output.tsv}

    # Find only primary ompA BLAST results (not secondary), use cat to avoid filename prefixes
    find {params.outdir} -path "*/best/blast/blast.ompa.tab" -type f ! -empty -exec cat {{}} \; >> {output.tsv}

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
    threads: 1
    shell: r"""
    echo -e "query\tsubject\tpident\tlength\tmismatch\tgapopen\tquery_start\tquery_end\tsubject_start\tsubject_end\tevalue\tbitscore" > {output.tsv}

    # Only add non-empty secondary BLAST files
    find {params.outdir} -path "*/best/blast/blast.secondary.tab" -type f ! -empty -exec cat {{}} \; >> {output.tsv}

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
