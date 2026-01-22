# workflow/rules/11-mixed-detection.smk
# Mixed strain detection via variant heterozygosity analysis

# Configuration
MIXED_DETECTION_ENABLED = config.get("mixed_detection", {}).get("enabled", True)
HET_THRESHOLD = config.get("mixed_detection", {}).get("het_threshold", 50)  # Sites to flag as mixed
MIN_AF = config.get("mixed_detection", {}).get("min_af", 0.15)  # Min allele freq for het site
MAX_AF = config.get("mixed_detection", {}).get("max_af", 0.85)  # Max allele freq for het site
MIN_DEPTH = config.get("mixed_detection", {}).get("min_depth", 10)  # Min read depth at site
MIN_BASE_QUAL = config.get("mixed_detection", {}).get("min_base_qual", 20)  # Min base quality
MIN_STRAND_READS = config.get("mixed_detection", {}).get("min_strand_reads", 2)  # Min reads on each strand

# -----------------------------------------------------------------------------
# INDEX BEST ASSEMBLY
# -----------------------------------------------------------------------------

rule mixed_extract_main_contig:
    """Extract the largest contig (N50) from assembly for mixed strain analysis.

    Focuses analysis on the main chromosome, reducing noise from small fragments,
    plasmid sequences, or assembly artifacts.
    """
    input:
        assembly = OUTDIR / "{sample}" / "best" / "assembly.fasta"
    output:
        main_contig = OUTDIR / "{sample}" / "mixed_detection" / "main_contig.fasta"
    conda: "../envs/bowtie.yaml"
    log: OUTDIR / "{sample}" / "log" / "mixed_detection.extract.{sample}.log"
    shell: """
        mkdir -p $(dirname {output.main_contig})

        # Find the largest contig and extract it
        # Using awk to find the longest sequence, then samtools to extract it
        largest_contig=$(awk '/^>/ {{if (seq) print length(seq), name; name=$1; seq=""}}
                             !/^>/ {{seq=seq$0}}
                             END {{print length(seq), name}}' {input.assembly} | \
                         sort -rn | head -1 | cut -d' ' -f2 | sed 's/>//')

        echo "Largest contig: $largest_contig" > {log}

        # Extract the largest contig using samtools
        samtools faidx {input.assembly}
        samtools faidx {input.assembly} "$largest_contig" > {output.main_contig} 2>> {log}
    """


rule mixed_index_assembly:
    """Index the main contig for read mapping."""
    input:
        main_contig = OUTDIR / "{sample}" / "mixed_detection" / "main_contig.fasta"
    output:
        index_done = OUTDIR / "{sample}" / "mixed_detection" / "index.done"
    params:
        prefix = lambda w: OUTDIR / w.sample / "mixed_detection" / "assembly"
    conda: "../envs/bowtie.yaml"
    log: OUTDIR / "{sample}" / "log" / "mixed_detection.index.{sample}.log"
    shell: """
        bowtie2-build {input.main_contig} {params.prefix} > {log} 2>&1
        touch {output.index_done}
    """

# -----------------------------------------------------------------------------
# MAP READS BACK TO ASSEMBLY
# -----------------------------------------------------------------------------

rule mixed_map_reads:
    """Map original reads back to the best assembly."""
    input:
        r1 = rules.extract_chlamydiales_reads.output.r1,
        r2 = rules.extract_chlamydiales_reads.output.r2,
        index_done = OUTDIR / "{sample}" / "mixed_detection" / "index.done"
    output:
        bam = OUTDIR / "{sample}" / "mixed_detection" / "{sample}.sorted.bam",
        bai = OUTDIR / "{sample}" / "mixed_detection" / "{sample}.sorted.bam.bai"
    params:
        index_prefix = lambda w: OUTDIR / w.sample / "mixed_detection" / "assembly"
    threads: config["threads"]["bowtie"]
    conda: "../envs/bowtie.yaml"
    log: OUTDIR / "{sample}" / "log" / "mixed_detection.map.{sample}.log"
    shell: """
        bowtie2 \
            -x {params.index_prefix} \
            -1 {input.r1} \
            -2 {input.r2} \
            --threads {threads} \
            2>> {log} | \
        samtools view -bSh - | \
        samtools sort -@ {threads} -o {output.bam} 2>> {log}

        samtools index {output.bam}
    """

# -----------------------------------------------------------------------------
# CALL VARIANTS AND ANALYZE ALLELE FREQUENCIES
# -----------------------------------------------------------------------------

rule mixed_call_variants:
    """Call variants and identify heterozygous-like sites indicating mixed strains.

    Analyzes only the main chromosome (largest contig) to reduce noise.

    Filters applied:
    - Allele frequency: 15-85% (intermediate, not fixed)
    - Minimum depth: 10x (configurable)
    - Base quality: ≥20 (reduces sequencing errors)
    - Strand bias: Requires reads on both strands (reduces systematic errors)
    """
    input:
        bam = OUTDIR / "{sample}" / "mixed_detection" / "{sample}.sorted.bam",
        main_contig = OUTDIR / "{sample}" / "mixed_detection" / "main_contig.fasta"
    output:
        vcf = OUTDIR / "{sample}" / "mixed_detection" / "{sample}.variants.vcf",
        het_sites = OUTDIR / "{sample}" / "mixed_detection" / "{sample}.het_sites.tsv",
        summary = OUTDIR / "{sample}" / "mixed_detection" / "{sample}.mixed_summary.txt",
        status = OUTDIR / "status" / "mixed_detection.{sample}.txt"
    params:
        min_af = MIN_AF,
        max_af = MAX_AF,
        min_depth = MIN_DEPTH,
        min_base_qual = MIN_BASE_QUAL,
        min_strand_reads = MIN_STRAND_READS,
        het_threshold = HET_THRESHOLD,
        sample = lambda w: w.sample
    threads: config["threads"]["samtools"]
    conda: "../envs/bowtie.yaml"
    log: OUTDIR / "{sample}" / "log" / "mixed_detection.variants.{sample}.log"
    shell: r"""
        # Index the main contig for bcftools
        samtools faidx {input.main_contig}

        # Call variants with bcftools (on main chromosome only)
        # -Q: minimum base quality filter
        # -a AD,ADF,ADR,DP: allele depth total, forward, reverse, and total depth
        bcftools mpileup \
            -f {input.main_contig} \
            -a AD,ADF,ADR,DP \
            -Q {params.min_base_qual} \
            --max-depth 1000 \
            {input.bam} 2>> {log} | \
        bcftools call \
            -mv \
            -Ov \
            -o {output.vcf} 2>> {log}

        # Extract heterozygous-like sites with strand bias filter
        echo -e "contig\tposition\tref\talt\tdepth\tref_count\talt_count\talt_freq\tfwd_alt\trev_alt\tstrand_bias" > {output.het_sites}

        bcftools query \
            -f '%CHROM\t%POS\t%REF\t%ALT\t%DP\t[%AD]\t[%ADF]\t[%ADR]\n' \
            {output.vcf} 2>> {log} | \
        awk -F'\t' -v min_af={params.min_af} -v max_af={params.max_af} -v min_dp={params.min_depth} -v min_strand={params.min_strand_reads} '
        BEGIN {{ OFS="\t" }}
        {{
            # Parse AD field (ref,alt counts)
            split($6, ad, ",")
            ref_count = ad[1]
            alt_count = ad[2]
            depth = $5

            # Parse strand-specific counts (ADF = forward, ADR = reverse)
            split($7, adf, ",")
            split($8, adr, ",")
            fwd_alt = (length(adf) > 1) ? adf[2] : 0
            rev_alt = (length(adr) > 1) ? adr[2] : 0

            if (depth >= min_dp && alt_count > 0) {{
                af = alt_count / depth

                # Check allele frequency range
                if (af >= min_af && af <= max_af) {{
                    # Check strand bias: require minimum reads on both strands
                    if (fwd_alt >= min_strand && rev_alt >= min_strand) {{
                        strand_bias = "PASS"
                        print $1, $2, $3, $4, depth, ref_count, alt_count, af, fwd_alt, rev_alt, strand_bias
                    }}
                }}
            }}
        }}' >> {output.het_sites}

        # Count heterozygous sites that passed all filters (excluding header)
        het_count=$(tail -n +2 {output.het_sites} | wc -l | tr -d ' ')

        # Calculate mean AF of het sites
        mean_af=$(tail -n +2 {output.het_sites} | awk '{{ sum += $8; n++ }} END {{ if (n > 0) printf "%.3f", sum/n; else print "N/A" }}')

        # Get the main contig name and size
        main_contig_name=$(head -1 {input.main_contig} | sed 's/>//' | cut -d' ' -f1)
        main_contig_size=$(tail -n +2 {input.main_contig} | tr -d '\n' | wc -c | tr -d ' ')

        # Check if variants are spread across chromosome or clustered
        if [ "$het_count" -gt 0 ]; then
            # Get range of positions across the main chromosome
            position_range=$(tail -n +2 {output.het_sites} | awk '
                BEGIN {{ min=999999999; max=0 }}
                {{ if ($2 < min) min=$2; if ($2 > max) max=$2 }}
                END {{ if (max > 0) print max - min; else print 0 }}
            ')
        else
            position_range=0
        fi

        # Generate summary report
        cat > {output.summary} << EOF
MIXED STRAIN DETECTION SUMMARY
==============================
Sample: {params.sample}
Analysis: Variant heterozygosity on main chromosome (largest contig)

PARAMETERS
----------
Min allele frequency: {params.min_af}
Max allele frequency: {params.max_af}
Min read depth: {params.min_depth}
Min base quality: {params.min_base_qual}
Min reads per strand: {params.min_strand_reads}
Het site threshold: {params.het_threshold}

RESULTS
-------
Main contig analyzed: $main_contig_name ($main_contig_size bp)
Heterozygous-like sites: $het_count
Mean allele frequency: $mean_af
Position range: $position_range bp

QUALITY FILTERS APPLIED
-----------------------
✓ Main chromosome only (excludes plasmid/fragments)
✓ Base quality ≥{params.min_base_qual} (reduces sequencing errors)
✓ Strand bias filter: ≥{params.min_strand_reads} reads on each strand (reduces systematic errors)
✓ Allele frequency {params.min_af}-{params.max_af} (excludes fixed variants and rare errors)

INTERPRETATION
--------------
EOF

        if [ "$het_count" -gt {params.het_threshold} ]; then
            echo "STATUS: MIXED_STRAIN_WARNING" >> {output.summary}
            echo "" >> {output.summary}
            echo "⚠️  HIGH HETEROZYGOSITY DETECTED" >> {output.summary}
            echo "Found $het_count sites with intermediate allele frequencies on main chromosome." >> {output.summary}
            echo "" >> {output.summary}
            if [ "$position_range" -gt 100000 ]; then
                echo "Distribution: GENOME-WIDE (variants span $position_range bp)" >> {output.summary}
                echo "Pattern consistent with true mixed infection." >> {output.summary}
            else
                echo "Distribution: LOCALIZED (variants within $position_range bp region)" >> {output.summary}
                echo "Could indicate recombination, hypervariable region, or artifact." >> {output.summary}
            fi
            echo "" >> {output.summary}
            echo "This suggests possible co-infection with multiple Ct strains." >> {output.summary}
            echo "" >> {output.summary}
            echo "Recommendations:" >> {output.summary}
            echo "  - Interpret typing results with caution" >> {output.summary}
            echo "  - ompA genovar may represent dominant strain only" >> {output.summary}
            echo "  - Consider strain-resolved analysis if clinically relevant" >> {output.summary}
        elif [ "$het_count" -gt 10 ]; then
            echo "STATUS: LOW_HETEROZYGOSITY" >> {output.summary}
            echo "" >> {output.summary}
            echo "Moderate heterozygosity detected ($het_count sites)." >> {output.summary}
            echo "Could indicate:" >> {output.summary}
            echo "  - Minor mixed population" >> {output.summary}
            echo "  - Within-host evolution" >> {output.summary}
            echo "  - Residual noise at problematic regions" >> {output.summary}
        else
            echo "STATUS: SINGLE_STRAIN" >> {output.summary}
            echo "" >> {output.summary}
            echo "Low heterozygosity ($het_count sites) consistent with single strain infection." >> {output.summary}
            echo "Quality filters effectively removed potential artifacts." >> {output.summary}
        fi

        touch {output.status}
    """

# -----------------------------------------------------------------------------
# COLLATE MIXED DETECTION RESULTS
# -----------------------------------------------------------------------------

rule mixed_detection_summary:
    """Collate mixed detection results across all samples."""
    input:
        summaries = expand(OUTDIR / "{sample}" / "mixed_detection" / "{sample}.mixed_summary.txt", sample=SAMPLES),
        het_sites = expand(OUTDIR / "{sample}" / "mixed_detection" / "{sample}.het_sites.tsv", sample=SAMPLES)
    output:
        summary_tsv = OUTDIR / "reports" / "mixed_detection_summary.tsv",
        report = OUTDIR / "reports" / "mixed_detection_report.txt",
        status = OUTDIR / "status" / "mixed_detection.summary.txt"
    params:
        samples = SAMPLES,
        outdir = OUTDIR,
        het_threshold = HET_THRESHOLD
    run:
        from pathlib import Path
        from datetime import datetime
        import re

        rows = []
        flagged_samples = []

        for sample in params.samples:
            summary_file = Path(params.outdir) / sample / "mixed_detection" / f"{sample}.mixed_summary.txt"
            het_file = Path(params.outdir) / sample / "mixed_detection" / f"{sample}.het_sites.tsv"

            row = {"sample": sample, "het_sites": 0, "mean_af": "N/A", "status": "Unknown"}

            if summary_file.exists():
                content = summary_file.read_text()

                # Extract het count
                match = re.search(r"Heterozygous-like sites:\s*(\d+)", content)
                if match:
                    row["het_sites"] = int(match.group(1))

                # Extract mean AF
                match = re.search(r"Mean allele frequency:\s*([\d.]+|N/A)", content)
                if match:
                    row["mean_af"] = match.group(1)

                # Extract status
                if "STATUS: MIXED_STRAIN_WARNING" in content:
                    row["status"] = "MIXED_WARNING"
                    flagged_samples.append(sample)
                elif "STATUS: LOW_HETEROZYGOSITY" in content:
                    row["status"] = "LOW_HET"
                elif "STATUS: SINGLE_STRAIN" in content:
                    row["status"] = "SINGLE"

            rows.append(row)

        # Write TSV summary
        with open(output.summary_tsv, 'w') as f:
            f.write("sample\thet_sites\tmean_af\tstatus\n")
            for row in sorted(rows, key=lambda x: x["sample"]):
                f.write(f"{row['sample']}\t{row['het_sites']}\t{row['mean_af']}\t{row['status']}\n")

        # Write text report
        with open(output.report, 'w') as f:
            f.write("=" * 70 + "\n")
            f.write("MIXED STRAIN DETECTION REPORT\n")
            f.write("=" * 70 + "\n")
            f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Total samples analyzed: {len(rows)}\n")
            f.write(f"Heterozygosity threshold: {params.het_threshold} sites\n\n")

            # Summary counts
            mixed_count = sum(1 for r in rows if r["status"] == "MIXED_WARNING")
            low_het_count = sum(1 for r in rows if r["status"] == "LOW_HET")
            single_count = sum(1 for r in rows if r["status"] == "SINGLE")

            f.write("-" * 40 + "\n")
            f.write("SUMMARY\n")
            f.write("-" * 40 + "\n")
            f.write(f"Single strain (low het):     {single_count} samples\n")
            f.write(f"Moderate heterozygosity:     {low_het_count} samples\n")
            f.write(f"⚠️  Mixed strain warning:     {mixed_count} samples\n\n")

            if flagged_samples:
                f.write("-" * 40 + "\n")
                f.write("SAMPLES WITH MIXED STRAIN WARNING\n")
                f.write("-" * 40 + "\n")
                for sample in sorted(flagged_samples):
                    row = next(r for r in rows if r["sample"] == sample)
                    f.write(f"  {sample}: {row['het_sites']} het sites (mean AF: {row['mean_af']})\n")
                f.write("\n")

            f.write("-" * 40 + "\n")
            f.write("ALL SAMPLES\n")
            f.write("-" * 40 + "\n")
            f.write(f"{'Sample':<30} {'Het Sites':>10} {'Mean AF':>10} {'Status':>15}\n")
            f.write("-" * 70 + "\n")
            for row in sorted(rows, key=lambda x: -x["het_sites"]):
                f.write(f"{row['sample']:<30} {row['het_sites']:>10} {row['mean_af']:>10} {row['status']:>15}\n")

            f.write("\n" + "=" * 70 + "\n")

        Path(output.status).touch()


# -----------------------------------------------------------------------------
# AGGREGATE RULE
# -----------------------------------------------------------------------------

rule all_mixed_detection:
    """Trigger mixed strain detection for all samples."""
    input:
        summary = OUTDIR / "reports" / "mixed_detection_summary.tsv",
        report = OUTDIR / "reports" / "mixed_detection_report.txt"
    output:
        status = OUTDIR / "status" / "mixed_detection.complete.txt"
    shell:
        "touch {output.status}"
