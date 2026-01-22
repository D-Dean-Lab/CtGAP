# workflow/rules/10-annotate.smk
# Genome annotation using Bakta for samples with sufficient coverage

import os
from pathlib import Path

# Configuration
ANNOTATION_ENABLED = config.get("annotation", {}).get("enabled", True)
MIN_COVERAGE_PCT = config.get("annotation", {}).get("min_coverage", 90)
CT_REFERENCE_SIZE = 1_043_000  # 1.043 Mb standard Ct genome size

# -----------------------------------------------------------------------------
# CHECKPOINT: Identify samples eligible for annotation
# -----------------------------------------------------------------------------

checkpoint identify_samples_for_annotation:
    """
    Identify samples with genome coverage >= threshold for annotation.

    Reads QUAST reports to determine assembly size, calculates coverage
    relative to standard Ct genome (1.04 Mb), and outputs list of passing samples.
    """
    input:
        quast_reports = expand(OUTDIR / "{sample}" / "denovo" / "assembly_statistics" / "report.tsv", sample=SAMPLES),
        best_assemblies = expand(OUTDIR / "{sample}" / "best" / "assembly.fasta", sample=SAMPLES)
    output:
        passing_samples = OUTDIR / "annotation" / "samples_to_annotate.txt",
        coverage_summary = OUTDIR / "annotation" / "coverage_check.tsv"
    params:
        samples = SAMPLES,
        outdir = OUTDIR,
        min_coverage = MIN_COVERAGE_PCT,
        ref_size = CT_REFERENCE_SIZE
    run:
        passing = []
        coverage_data = []

        for sample in params.samples:
            quast_file = Path(params.outdir) / sample / "denovo" / "assembly_statistics" / "report.tsv"
            n50 = 0

            if quast_file.exists():
                with open(quast_file) as f:
                    for line in f:
                        # Use N50 as it represents contiguous assembly after scaffolding
                        if line.startswith("N50"):
                            parts = line.strip().split("\t")
                            if len(parts) >= 2:
                                n50 = int(parts[1].replace(",", ""))
                            break

            coverage_pct = (n50 / params.ref_size) * 100 if n50 > 0 else 0
            passes = coverage_pct >= params.min_coverage

            coverage_data.append({
                "sample": sample,
                "n50": n50,
                "coverage_pct": f"{coverage_pct:.1f}",
                "passes_threshold": "YES" if passes else "NO"
            })

            if passes:
                passing.append(sample)

        # Write passing samples list
        with open(output.passing_samples, 'w') as f:
            f.write("\n".join(sorted(passing)))

        # Write coverage summary TSV
        with open(output.coverage_summary, 'w') as f:
            f.write("sample\tn50\tcoverage_pct\tpasses_threshold\n")
            for row in sorted(coverage_data, key=lambda x: x["sample"]):
                f.write(f"{row['sample']}\t{row['n50']}\t{row['coverage_pct']}\t{row['passes_threshold']}\n")

        # Log summary
        print(f"Annotation checkpoint: {len(passing)}/{len(list(params.samples))} samples pass {params.min_coverage}% coverage threshold")


# -----------------------------------------------------------------------------
# INPUT FUNCTION: Get samples that pass annotation threshold
# -----------------------------------------------------------------------------

def get_annotation_inputs(wildcards):
    """
    Input function that reads checkpoint output to determine which samples
    should be annotated. Returns list of Bakta output files for passing samples.
    """
    checkpoint_output = checkpoints.identify_samples_for_annotation.get(**wildcards).output.passing_samples

    with open(checkpoint_output) as f:
        samples = [s.strip() for s in f.readlines() if s.strip()]

    return expand(OUTDIR / "{sample}" / "annotation" / "{sample}.gff3", sample=samples)


def get_annotation_sample_list(wildcards):
    """
    Returns list of samples that passed annotation threshold.
    """
    checkpoint_output = checkpoints.identify_samples_for_annotation.get(**wildcards).output.passing_samples

    with open(checkpoint_output) as f:
        samples = [s.strip() for s in f.readlines() if s.strip()]

    return samples


# -----------------------------------------------------------------------------
# BAKTA ANNOTATION RULE
# -----------------------------------------------------------------------------

rule bakta_annotate:
    """
    Annotate genome assembly using Bakta.

    Only runs for samples that pass the coverage threshold (determined by checkpoint).
    Produces GFF3, GenBank, protein FASTA, and other standard annotation outputs.
    Requires: Bakta database downloaded via 'ctgap setup'.
    """
    input:
        assembly = OUTDIR / "{sample}" / "best" / "assembly.fasta"
    output:
        gff = OUTDIR / "{sample}" / "annotation" / "{sample}.gff3",
        gbk = OUTDIR / "{sample}" / "annotation" / "{sample}.gbff",
        faa = OUTDIR / "{sample}" / "annotation" / "{sample}.faa",
        ffn = OUTDIR / "{sample}" / "annotation" / "{sample}.ffn",
        tsv = OUTDIR / "{sample}" / "annotation" / "{sample}.tsv",
        txt = OUTDIR / "{sample}" / "annotation" / "{sample}.txt",
        status = OUTDIR / "status" / "annotation.{sample}.txt"
    params:
        db = config.get("bakta", {}).get("db", "resources/bakta_db/db-light"),
        genus = config.get("bakta", {}).get("genus", "Chlamydia"),
        species = config.get("bakta", {}).get("species", "trachomatis"),
        gram = config.get("bakta", {}).get("gram", "?"),
        outdir = lambda w: OUTDIR / w.sample / "annotation",
        prefix = lambda w: w.sample
    threads: config["threads"]["bakta"]
    conda: "../envs/bakta.yaml"
    log: OUTDIR / "{sample}" / "log" / "annotation.{sample}.log"
    benchmark: OUTDIR / "{sample}" / "benchmark" / "annotation.{sample}.txt"
    shell: """
        bakta \
            --db {params.db} \
            --output {params.outdir} \
            --prefix {params.prefix} \
            --genus {params.genus} \
            --species {params.species} \
            --gram {params.gram} \
            --threads {threads} \
            --force \
            {input.assembly} \
            > {log} 2>&1

        touch {output.status}
    """


# -----------------------------------------------------------------------------
# ANNOTATION SUMMARY
# -----------------------------------------------------------------------------

rule annotation_summary:
    """
    Create summary of annotation results across all annotated samples.
    """
    input:
        annotations = get_annotation_inputs,
        coverage_summary = OUTDIR / "annotation" / "coverage_check.tsv",
        passing_samples = OUTDIR / "annotation" / "samples_to_annotate.txt"
    output:
        summary = OUTDIR / "reports" / "annotation_summary.tsv",
        report = OUTDIR / "reports" / "annotation_report.txt",
        status = OUTDIR / "status" / "annotation.summary.txt"
    params:
        outdir = OUTDIR
    run:
        import re
        from datetime import datetime

        # Read samples from checkpoint output file
        with open(input.passing_samples) as f:
            samples = [s.strip() for s in f.readlines() if s.strip()]
        rows = []

        for sample in samples:
            txt_file = Path(params.outdir) / sample / "annotation" / f"{sample}.txt"
            row = {"sample": sample}

            if txt_file.exists():
                content = txt_file.read_text()

                # Parse Bakta summary statistics
                # Look for patterns like "CDSs: 920" or "tRNAs: 37"
                patterns = {
                    "contigs": r"Sequences:\s*(\d+)",
                    "size_bp": r"Size:\s*([\d,]+)",
                    "gc_percent": r"GC:\s*([\d.]+)",
                    "cds_count": r"CDSs:\s*(\d+)",
                    "trna_count": r"tRNAs:\s*(\d+)",
                    "rrna_count": r"rRNAs:\s*(\d+)",
                    "ncrna_count": r"ncRNAs:\s*(\d+)",
                    "crispr_count": r"CRISPRs:\s*(\d+)",
                    "hypothetical_count": r"hypotheticals:\s*(\d+)",
                    "pseudogene_count": r"pseudogenes:\s*(\d+)"
                }

                for key, pattern in patterns.items():
                    match = re.search(pattern, content, re.IGNORECASE)
                    if match:
                        val = match.group(1).replace(",", "")
                        row[key] = val
                    else:
                        row[key] = ""

            rows.append(row)

        # Write TSV summary
        fieldnames = ["sample", "contigs", "size_bp", "gc_percent", "cds_count",
                      "trna_count", "rrna_count", "ncrna_count", "hypothetical_count", "pseudogene_count"]

        with open(output.summary, "w") as f:
            f.write("\t".join(fieldnames) + "\n")
            for row in sorted(rows, key=lambda x: x["sample"]):
                f.write("\t".join(str(row.get(k, "")) for k in fieldnames) + "\n")

        # Write text report
        with open(output.report, "w") as f:
            f.write("=" * 70 + "\n")
            f.write("CtGAP GENOME ANNOTATION SUMMARY (Bakta)\n")
            f.write("=" * 70 + "\n")
            f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Samples annotated: {len(rows)}\n")
            f.write(f"Minimum coverage threshold: {MIN_COVERAGE_PCT}%\n\n")

            # Read coverage summary
            f.write("-" * 40 + "\n")
            f.write("COVERAGE CHECK RESULTS\n")
            f.write("-" * 40 + "\n")
            with open(input.coverage_summary) as cs:
                f.write(cs.read())
            f.write("\n")

            if rows:
                # Summary statistics
                cds_counts = [int(r["cds_count"]) for r in rows if r.get("cds_count", "").isdigit()]

                f.write("-" * 40 + "\n")
                f.write("ANNOTATION STATISTICS\n")
                f.write("-" * 40 + "\n")

                if cds_counts:
                    f.write(f"CDS count range: {min(cds_counts)} - {max(cds_counts)}\n")
                    f.write(f"Median CDS count: {sorted(cds_counts)[len(cds_counts)//2]}\n")

                f.write("\n")
                f.write("-" * 40 + "\n")
                f.write("PER-SAMPLE DETAILS\n")
                f.write("-" * 40 + "\n")

                for row in sorted(rows, key=lambda x: x["sample"]):
                    f.write(f"\n{row['sample']}:\n")
                    f.write(f"  CDSs: {row.get('cds_count', 'N/A')}\n")
                    f.write(f"  tRNAs: {row.get('trna_count', 'N/A')}\n")
                    f.write(f"  rRNAs: {row.get('rrna_count', 'N/A')}\n")
                    f.write(f"  Hypotheticals: {row.get('hypothetical_count', 'N/A')}\n")
            else:
                f.write("No samples met the coverage threshold for annotation.\n")

            f.write("\n" + "=" * 70 + "\n")

        Path(output.status).touch()


# -----------------------------------------------------------------------------
# AGGREGATE RULE (triggers annotation pipeline)
# -----------------------------------------------------------------------------

rule all_annotations:
    """
    Aggregate rule that triggers annotation for all eligible samples.
    Used as a target in the main workflow.
    """
    input:
        summary = OUTDIR / "reports" / "annotation_summary.tsv",
        report = OUTDIR / "reports" / "annotation_report.txt"
    output:
        status = OUTDIR / "status" / "annotation.complete.txt"
    shell:
        "touch {output.status}"
