# =============================================================================
# DEDICATED PLASMID ASSEMBLY MODULE
# =============================================================================
# Reference-guided read extraction + Unicycler assembly for C. trachomatis plasmid
#
# Strategy:
#   1. Map reads to reference plasmid (D strain pCT, ~7,500 bp)
#   2. Extract mapped read pairs (plasmid-specific reads)
#   3. Assemble extracted reads with Unicycler
#   4. BLAST verify and type the assembly
#
# This approach is cleaner than filtering whole-genome assemblies by size,
# as it focuses the assembler on plasmid reads only.
#
# C. trachomatis plasmid facts:
#   - Highly conserved (>95% identity across strains)
#   - ~7,500 bp circular
#   - Present in most CT strains (some are plasmid-free)
# =============================================================================

import os
from pathlib import Path

# Reference plasmid (D strain - standard CT reference)
REF_PLASMID = RESOURCES / "references" / ORG / "reference_plasmid.fasta"

# Plasmid size expectations
CT_PLASMID_EXPECTED = 7500
CT_PLASMID_MIN = 6000
CT_PLASMID_MAX = 10000

# -----------------------------------------------------------------------------
# 1. INDEX REFERENCE PLASMID
# -----------------------------------------------------------------------------
rule index_reference_plasmid:
    """Create bowtie2 index for reference plasmid."""
    input:
        ref = REF_PLASMID
    output:
        idx = multiext(str(OUTDIR / "plasmid" / "ref_index" / "pCT"),
                      ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"),
        done = OUTDIR / "status" / "plasmid.ref_index.txt"
    params:
        prefix = OUTDIR / "plasmid" / "ref_index" / "pCT"
    conda: "../envs/bowtie.yaml"
    log: OUTDIR / "plasmid" / "log" / "index_reference_plasmid.log"
    threads: config["threads"]["bowtieindex"]
    shell: r"""
        mkdir -p $(dirname {params.prefix})
        bowtie2-build {input.ref} {params.prefix} > {log} 2>&1
        touch {output.done}
    """

# -----------------------------------------------------------------------------
# 2. MAP READS TO REFERENCE PLASMID
# -----------------------------------------------------------------------------
rule map_reads_to_plasmid:
    """Map scrubbed reads to reference plasmid to identify plasmid reads."""
    input:
        r1 = OUTDIR / "{sample}" / "scrub" / "{sample}_trim_scrub_R1.fastq.gz",
        r2 = OUTDIR / "{sample}" / "scrub" / "{sample}_trim_scrub_R2.fastq.gz",
        idx_done = OUTDIR / "status" / "plasmid.ref_index.txt"
    output:
        bam = OUTDIR / "{sample}" / "plasmid" / "mapping" / "{sample}.plasmid_mapped.bam",
        stats = OUTDIR / "{sample}" / "plasmid" / "mapping" / "{sample}.mapping_stats.txt"
    params:
        idx_prefix = OUTDIR / "plasmid" / "ref_index" / "pCT"
    conda: "../envs/bowtie.yaml"
    log: OUTDIR / "{sample}" / "log" / "plasmid.map_reads.{sample}.log"
    threads: config["threads"]["bowtie"]
    shell: r"""
        mkdir -p $(dirname {output.bam})

        # Map reads to reference plasmid
        # Use sensitive-local for better capture of divergent regions
        bowtie2 \
            -x {params.idx_prefix} \
            -1 {input.r1} \
            -2 {input.r2} \
            --threads {threads} \
            --sensitive-local \
            2> {log} | \
        samtools view -bS -F 4 - | \
        samtools sort -@ {threads} -o {output.bam} -

        samtools index {output.bam}

        # Generate mapping stats
        echo "Plasmid Read Mapping Statistics" > {output.stats}
        echo "Sample: {wildcards.sample}" >> {output.stats}
        echo "Reference: pCT (D strain)" >> {output.stats}
        echo "================================" >> {output.stats}
        samtools flagstat {output.bam} >> {output.stats}
    """

# -----------------------------------------------------------------------------
# 3. EXTRACT PLASMID READ PAIRS
# -----------------------------------------------------------------------------
rule extract_plasmid_reads:
    """Extract read pairs that mapped to the reference plasmid."""
    input:
        bam = OUTDIR / "{sample}" / "plasmid" / "mapping" / "{sample}.plasmid_mapped.bam"
    output:
        r1 = OUTDIR / "{sample}" / "plasmid" / "reads" / "{sample}.plasmid_R1.fastq.gz",
        r2 = OUTDIR / "{sample}" / "plasmid" / "reads" / "{sample}.plasmid_R2.fastq.gz",
        read_count = OUTDIR / "{sample}" / "plasmid" / "reads" / "{sample}.read_count.txt"
    conda: "../envs/bowtie.yaml"
    log: OUTDIR / "{sample}" / "log" / "plasmid.extract_reads.{sample}.log"
    threads: config["threads"]["samtools"]
    shell: r"""
        mkdir -p $(dirname {output.r1})

        # Extract properly paired reads
        # -f 3: both reads mapped, proper pair
        samtools view -bf 3 {input.bam} | \
        samtools sort -n -@ {threads} - | \
        samtools fastq -@ {threads} \
            -1 {output.r1} \
            -2 {output.r2} \
            -0 /dev/null \
            -s /dev/null \
            -n \
            2> {log}

        # Count extracted reads (gzip -dc works on both Mac and Linux)
        num_reads=$(gzip -dc {output.r1} | wc -l | awk '{{print $1/4}}')
        echo "Extracted plasmid read pairs: $num_reads" > {output.read_count}
        echo "$num_reads" >> {output.read_count}
    """

# -----------------------------------------------------------------------------
# 4. UNICYCLER PLASMID ASSEMBLY
# -----------------------------------------------------------------------------
rule unicycler_plasmid_assembly:
    """Assemble extracted plasmid reads with Unicycler."""
    input:
        r1 = OUTDIR / "{sample}" / "plasmid" / "reads" / "{sample}.plasmid_R1.fastq.gz",
        r2 = OUTDIR / "{sample}" / "plasmid" / "reads" / "{sample}.plasmid_R2.fastq.gz",
        count_file = OUTDIR / "{sample}" / "plasmid" / "reads" / "{sample}.read_count.txt"
    output:
        assembly = OUTDIR / "{sample}" / "plasmid" / "assembly" / "assembly.fasta",
        graph = OUTDIR / "{sample}" / "plasmid" / "assembly" / "assembly.gfa",
        log_file = OUTDIR / "{sample}" / "plasmid" / "assembly" / "unicycler.log"
    params:
        outdir = lambda w: str(OUTDIR / w.sample / "plasmid" / "assembly"),
        min_reads = 100  # Minimum reads needed for assembly
    conda: "../envs/unicycler.yaml"
    log: OUTDIR / "{sample}" / "log" / "plasmid.unicycler.{sample}.log"
    threads: config['threads']['spades']
    shell: r"""
        mkdir -p {params.outdir}

        # Check if we have enough reads
        read_count=$(tail -1 {input.count_file})

        if [ "$read_count" -lt "{params.min_reads}" ]; then
            echo "Insufficient plasmid reads ($read_count < {params.min_reads})" > {log}
            echo "Sample may be plasmid-free or have low plasmid coverage" >> {log}
            echo ">no_plasmid_detected" > {output.assembly}
            echo "NNNNNNNNNN" >> {output.assembly}
            touch {output.graph}
            echo "PLASMID_FREE_OR_LOW_COVERAGE" > {output.log_file}
        else
            # Run Unicycler
            unicycler \
                -1 {input.r1} \
                -2 {input.r2} \
                -o {params.outdir} \
                -t {threads} \
                --min_fasta_length 1000 \
                --verbosity 2 \
                > {log} 2>&1 || {{
                    echo "Unicycler failed" >> {log}
                    echo ">assembly_failed" > {output.assembly}
                    echo "NNNNNNNNNN" >> {output.assembly}
                    touch {output.graph}
                    echo "UNICYCLER_FAILED" > {output.log_file}
                }}
        fi
    """

# -----------------------------------------------------------------------------
# 5. ANALYZE PLASMID ASSEMBLY
# -----------------------------------------------------------------------------
rule analyze_plasmid_assembly:
    """Analyze the plasmid assembly - check size, circularity, completeness."""
    input:
        assembly = OUTDIR / "{sample}" / "plasmid" / "assembly" / "assembly.fasta",
        graph = OUTDIR / "{sample}" / "plasmid" / "assembly" / "assembly.gfa",
        mapping_stats = OUTDIR / "{sample}" / "plasmid" / "mapping" / "{sample}.mapping_stats.txt",
        read_count = OUTDIR / "{sample}" / "plasmid" / "reads" / "{sample}.read_count.txt"
    output:
        stats = OUTDIR / "{sample}" / "plasmid" / "{sample}.plasmid_analysis.txt",
        final_plasmid = OUTDIR / "{sample}" / "plasmid" / "{sample}.plasmid.fasta",
        status = OUTDIR / "status" / "plasmid.analysis.{sample}.txt"
    params:
        sample = lambda w: w.sample,
        expected_size = CT_PLASMID_EXPECTED,
        min_size = CT_PLASMID_MIN,
        max_size = CT_PLASMID_MAX
    shell: r"""
        python3 << 'PYTHON_SCRIPT'
import os
import sys

sample = "{params.sample}"
expected_size = {params.expected_size}
min_size = {params.min_size}
max_size = {params.max_size}

assembly_file = "{input.assembly}"
graph_file = "{input.graph}"
stats_file = "{output.stats}"
final_plasmid = "{output.final_plasmid}"

# Parse assembly
contigs = []
with open(assembly_file) as f:
    header = ""
    seq = ""
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            if header and seq:
                contigs.append((header, seq, len(seq)))
            header = line[1:]
            seq = ""
        else:
            seq += line
    if header and seq:
        contigs.append((header, seq, len(seq)))

# Check for circular contigs in GFA
circular_contigs = set()
if os.path.exists(graph_file) and os.path.getsize(graph_file) > 0:
    with open(graph_file) as f:
        for line in f:
            if line.startswith("L") and "circular" in line.lower():
                circular_contigs.add(line.split()[1])

# Define newline character (avoids escaping issues in Snakemake heredoc)
NL = chr(10)

# Analyze and write stats
with open(stats_file, "w") as sf, open(final_plasmid, "w") as pf:
    sf.write("=" * 60 + NL)
    sf.write("PLASMID ANALYSIS REPORT" + NL)
    sf.write("=" * 60 + NL)
    sf.write("Sample: " + sample + NL)
    sf.write("Expected CT plasmid size: ~" + str(expected_size) + " bp" + NL)
    sf.write("Acceptable range: " + str(min_size) + "-" + str(max_size) + " bp" + NL)
    sf.write(NL)

    # Check for failed assembly
    if len(contigs) == 1 and contigs[0][0] in ["no_plasmid_detected", "assembly_failed"]:
        sf.write("STATUS: NO PLASMID DETECTED" + NL)
        sf.write(NL + "Possible reasons:" + NL)
        sf.write("  - Sample is plasmid-free (some CT strains lack plasmids)" + NL)
        sf.write("  - Low sequencing depth of plasmid" + NL)
        sf.write("  - Plasmid reads did not map to reference" + NL)
        pf.write(">no_plasmid" + NL + "NNNNNNNNNN" + NL)
    else:
        sf.write("CONTIGS FOUND: " + str(len(contigs)) + NL)
        sf.write("-" * 40 + NL)

        plasmid_found = False
        best_plasmid = None
        best_score = 0

        for header, seq, length in contigs:
            # Determine classification
            is_circular = header.split()[0] in circular_contigs or "circular" in header.lower()

            if min_size <= length <= max_size:
                classification = "CT_PLASMID"
                confidence = "HIGH"
                # Score: closer to expected size = better
                score = 100 - abs(length - expected_size) / expected_size * 100
                if is_circular:
                    score += 20  # Bonus for circular
                    confidence = "VERY HIGH"

                if score > best_score:
                    best_score = score
                    best_plasmid = (header, seq, length, is_circular)

                plasmid_found = True
            elif length < min_size:
                classification = "TOO_SMALL"
                confidence = "N/A"
            elif length > max_size:
                classification = "TOO_LARGE"
                confidence = "N/A"
            else:
                classification = "UNKNOWN"
                confidence = "N/A"

            circular_str = " [CIRCULAR]" if is_circular else ""
            sf.write("  " + header.split()[0] + ": " + str(length) + " bp" + circular_str + " - " + classification + NL)

        sf.write(NL)

        if best_plasmid:
            header, seq, length, is_circular = best_plasmid
            sf.write("RESULT: CT PLASMID DETECTED" + NL)
            sf.write("  Length: " + str(length) + " bp" + NL)
            sf.write("  Circular: " + ("Yes" if is_circular else "No/Unknown") + NL)
            sf.write("  Size deviation: " + str(abs(length - expected_size)) + " bp from expected" + NL)

            # Write final plasmid
            pf.write(">" + sample + "_plasmid length=" + str(length))
            if is_circular:
                pf.write(" topology=circular")
            pf.write(NL)
            # Write sequence in 80-char lines
            for i in range(0, len(seq), 80):
                pf.write(seq[i:i+80] + NL)
        else:
            sf.write("RESULT: NO VALID CT PLASMID FOUND" + NL)
            pf.write(">no_valid_plasmid" + NL + "NNNNNNNNNN" + NL)

print("Plasmid analysis complete for " + sample)
PYTHON_SCRIPT

        touch {output.status}
    """

# -----------------------------------------------------------------------------
# 5b. RE-ORIENT PLASMID TO MATCH REFERENCE
# -----------------------------------------------------------------------------
rule reorient_plasmid:
    """
    Re-orient plasmid to match reference orientation and start position.

    Circular plasmids may be assembled in reverse complement or with a different
    start position. This causes BLAST to report split alignments. This rule:
    1. Aligns plasmid to reference using minimap2
    2. Reverse complements if on minus strand
    3. Rotates sequence to align start positions
    """
    input:
        plasmid = OUTDIR / "{sample}" / "plasmid" / "{sample}.plasmid.fasta",
        reference = REF_PLASMID
    output:
        reoriented = OUTDIR / "{sample}" / "plasmid" / "{sample}.plasmid.reoriented.fasta",
        orientation_log = OUTDIR / "{sample}" / "plasmid" / "{sample}.orientation.txt"
    conda: "../envs/bowtie.yaml"
    log: OUTDIR / "{sample}" / "log" / "plasmid.reorient.{sample}.log"
    script: "../scripts/reorient_plasmid.py"


# -----------------------------------------------------------------------------
# 6. BLAST PLASMID AGAINST REFERENCES
# -----------------------------------------------------------------------------
rule blast_plasmid_dedicated:
    """BLAST the assembled plasmid against reference plasmids."""
    input:
        plasmid = OUTDIR / "{sample}" / "plasmid" / "{sample}.plasmid.reoriented.fasta",
        blastdb = "resources/references/ct/20_Ct_plasmids.nhr"
    output:
        blast = OUTDIR / "{sample}" / "plasmid" / "{sample}.plasmid_blast.tsv",
        status = OUTDIR / "status" / "plasmid.blast.{sample}.txt"
    conda: "../envs/misc.yaml"
    log: OUTDIR / "{sample}" / "log" / "plasmid.blast.{sample}.log"
    threads: config["threads"]["blast"]
    shell: r"""
        # Check if plasmid is valid (not placeholder)
        if grep -q "no_plasmid\|no_valid_plasmid\|NNNNNNNNNN" {input.plasmid}; then
            echo "No valid plasmid to BLAST" > {log}
            echo -e "query\tsubject\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore" > {output.blast}
        else
            blastn \
                -query {input.plasmid} \
                -db resources/references/ct/20_Ct_plasmids \
                -outfmt 6 \
                -max_target_seqs 5 \
                -out {output.blast} \
                2> {log}
        fi

        touch {output.status}
    """

# -----------------------------------------------------------------------------
# 7. MLST PLASMID TYPING
# -----------------------------------------------------------------------------
rule mlst_plasmid_dedicated:
    """MLST typing of the assembled plasmid."""
    input:
        plasmid = OUTDIR / "{sample}" / "plasmid" / "{sample}.plasmid.fasta"
    output:
        mlst = OUTDIR / "{sample}" / "plasmid" / "{sample}.plasmid_mlst.txt",
        status = OUTDIR / "status" / "plasmid.mlst.{sample}.txt"
    conda: "../envs/mlst.yaml"
    log: OUTDIR / "{sample}" / "log" / "plasmid.mlst.{sample}.log"
    shell: r"""
        # Check if plasmid is valid
        if grep -q "no_plasmid\|no_valid_plasmid\|NNNNNNNNNN" {input.plasmid}; then
            echo "No valid plasmid for MLST" > {log}
            echo -e "Sample\tST\tAlleles" > {output.mlst}
            echo -e "{wildcards.sample}\tNA\tNA" >> {output.mlst}
        else
            claMLST search \
                resources/references/ct/pubmlst/plasmid \
                {input.plasmid} > {output.mlst} 2> {log} || {{
                    echo -e "Sample\tST\tAlleles" > {output.mlst}
                    echo -e "{wildcards.sample}\tNA\tNA" >> {output.mlst}
                }}
        fi

        touch {output.status}
    """

# -----------------------------------------------------------------------------
# 8. COLLATE ALL PLASMID RESULTS
# -----------------------------------------------------------------------------
rule collate_plasmid_results:
    """Collate plasmid analysis results across all samples."""
    input:
        analyses = expand(OUTDIR / "{sample}" / "plasmid" / "{sample}.plasmid_analysis.txt", sample=SAMPLES),
        blasts = expand(OUTDIR / "{sample}" / "plasmid" / "{sample}.plasmid_blast.tsv", sample=SAMPLES),
        mlsts = expand(OUTDIR / "{sample}" / "plasmid" / "{sample}.plasmid_mlst.txt", sample=SAMPLES)
    output:
        summary = OUTDIR / "reports" / "plasmid_summary.tsv",
        blast_all = OUTDIR / "plasmid.blast.all_samples.tsv",
        mlst_all = OUTDIR / "plasmid.mlst.all_samples.tsv",
        report = OUTDIR / "reports" / "plasmid_analysis_report.txt",
        status = OUTDIR / "status" / "plasmid.collate.txt"
    params:
        samples = SAMPLES,
        outdir = OUTDIR
    run:
        import csv
        from pathlib import Path

        results = []

        for sample in params.samples:
            row = {"sample": sample}

            # Parse analysis file
            analysis_file = Path(params.outdir) / sample / "plasmid" / f"{sample}.plasmid_analysis.txt"
            if analysis_file.exists():
                content = analysis_file.read_text()
                if "CT PLASMID DETECTED" in content:
                    row["plasmid_status"] = "DETECTED"
                    # Extract length
                    for line in content.split("\n"):
                        if "Length:" in line:
                            row["plasmid_length"] = line.split(":")[1].strip().replace(" bp", "")
                        if "Circular:" in line:
                            row["circular"] = line.split(":")[1].strip()
                elif "NO PLASMID DETECTED" in content:
                    row["plasmid_status"] = "NOT_DETECTED"
                    row["plasmid_length"] = ""
                    row["circular"] = ""
                else:
                    row["plasmid_status"] = "NO_VALID_PLASMID"
                    row["plasmid_length"] = ""
                    row["circular"] = ""

            # Parse BLAST file for best hit
            blast_file = Path(params.outdir) / sample / "plasmid" / f"{sample}.plasmid_blast.tsv"
            row["blast_hit"] = ""
            row["blast_identity"] = ""
            if blast_file.exists() and blast_file.stat().st_size > 0:
                with open(blast_file) as f:
                    for line in f:
                        if not line.startswith("query"):
                            fields = line.strip().split("\t")
                            if len(fields) >= 3:
                                row["blast_hit"] = fields[1]
                                row["blast_identity"] = fields[2]
                            break

            results.append(row)

        # Write summary TSV
        with open(output.summary, "w", newline="") as f:
            fieldnames = ["sample", "plasmid_status", "plasmid_length", "circular", "blast_hit", "blast_identity"]
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
            writer.writeheader()
            writer.writerows(results)

        # Collate all BLAST results
        with open(output.blast_all, "w") as f:
            f.write("query\tsubject\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\n")
            for sample in params.samples:
                blast_file = Path(params.outdir) / sample / "plasmid" / f"{sample}.plasmid_blast.tsv"
                if blast_file.exists():
                    with open(blast_file) as bf:
                        for line in bf:
                            if not line.startswith("query"):
                                f.write(line)

        # Collate all MLST results
        with open(output.mlst_all, "w") as f:
            first = True
            for sample in params.samples:
                mlst_file = Path(params.outdir) / sample / "plasmid" / f"{sample}.plasmid_mlst.txt"
                if mlst_file.exists():
                    with open(mlst_file) as mf:
                        for i, line in enumerate(mf):
                            if first or i > 0:  # Include header only once
                                f.write(line)
                        first = False

        # Generate report
        with open(output.report, "w") as f:
            detected = sum(1 for r in results if r.get("plasmid_status") == "DETECTED")
            not_detected = sum(1 for r in results if r.get("plasmid_status") == "NOT_DETECTED")

            f.write("=" * 60 + "\n")
            f.write("CT PLASMID ANALYSIS SUMMARY\n")
            f.write("=" * 60 + "\n\n")
            f.write(f"Total samples: {len(results)}\n")
            f.write(f"Plasmid detected: {detected}\n")
            f.write(f"Plasmid not detected: {not_detected}\n")
            f.write(f"Other: {len(results) - detected - not_detected}\n\n")

            f.write("-" * 40 + "\n")
            f.write("METHODOLOGY\n")
            f.write("-" * 40 + "\n")
            f.write("1. Reads mapped to reference plasmid (pCT, D strain)\n")
            f.write("2. Plasmid-specific read pairs extracted\n")
            f.write("3. Unicycler assembly of extracted reads\n")
            f.write("4. Size verification (expected ~7,500 bp)\n")
            f.write("5. BLAST confirmation against 20 CT plasmids\n\n")

            f.write("-" * 40 + "\n")
            f.write("PER-SAMPLE RESULTS\n")
            f.write("-" * 40 + "\n")
            for r in results:
                f.write(f"\n{r['sample']}:\n")
                f.write(f"  Status: {r.get('plasmid_status', 'UNKNOWN')}\n")
                if r.get('plasmid_length'):
                    f.write(f"  Length: {r['plasmid_length']} bp\n")
                if r.get('circular'):
                    f.write(f"  Circular: {r['circular']}\n")
                if r.get('blast_hit'):
                    f.write(f"  Best BLAST hit: {r['blast_hit']} ({r['blast_identity']}% identity)\n")

        Path(output.status).touch()

# -----------------------------------------------------------------------------
# 9. PLASMID PDF REPORT (per sample)
# -----------------------------------------------------------------------------
rule plasmid_pdf_report_dedicated:
    """Generate PDF report for plasmid analysis."""
    input:
        analysis = OUTDIR / "{sample}" / "plasmid" / "{sample}.plasmid_analysis.txt",
        plasmid = OUTDIR / "{sample}" / "plasmid" / "{sample}.plasmid.fasta",
        blast = OUTDIR / "{sample}" / "plasmid" / "{sample}.plasmid_blast.tsv",
        mlst = OUTDIR / "{sample}" / "plasmid" / "{sample}.plasmid_mlst.txt",
        mapping_stats = OUTDIR / "{sample}" / "plasmid" / "mapping" / "{sample}.mapping_stats.txt"
    output:
        pdf = OUTDIR / "{sample}" / "reports" / "{sample}_plasmid_report.pdf"
    conda: "../envs/pdf.yaml"
    log: OUTDIR / "{sample}" / "log" / "plasmid.pdf_report.{sample}.log"
    script: "../scripts/create_plasmid_report.py"
