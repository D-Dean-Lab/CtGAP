# workflow/rules/6b-filter-contigs.smk
# Post-assembly filtering to remove non-Ct contigs

# Configuration
FILTER_ENABLED = config.get("assembly_filter", {}).get("enabled", True)
MIN_IDENTITY = config.get("assembly_filter", {}).get("min_identity", 90)
MIN_COVERAGE = config.get("assembly_filter", {}).get("min_coverage", 50)

# Reference for filtering (use plurality consensus - tolerates strain variation)
FILTER_REF = RESOURCES / "references" / ORG / "plurality.fasta"

if FILTER_ENABLED:

    rule filter_assembly_contigs:
        """
        Filter assembly to remove non-Ct contigs.

        Uses minimap2 to align contigs to Ct reference and keeps only those
        with sufficient identity and coverage. This removes contamination
        and assembly artifacts.
        """
        input:
            assembly = OUTDIR / "{sample}" / "best" / "assembly.fasta",
            method = OUTDIR / "reports" / "assembly_selection" / "{sample}.method.txt"
        output:
            filtered = OUTDIR / "{sample}" / "best" / "assembly.filtered.fasta",
            report = OUTDIR / "{sample}" / "best" / "filter_report.txt",
            status = OUTDIR / "status" / "filter.{sample}.txt"
        params:
            reference = FILTER_REF,
            min_identity = MIN_IDENTITY,
            min_coverage = MIN_COVERAGE,
            sample = lambda w: w.sample
        conda: "../envs/bowtie.yaml"
        log: OUTDIR / "{sample}" / "log" / "filter_contigs.{sample}.log"
        shell: r"""
            set -e

            # Create output directory
            mkdir -p $(dirname {output.filtered})

            # Get contig names and lengths from assembly
            samtools faidx {input.assembly}

            # Align contigs to Ct reference
            minimap2 -x asm5 -t 4 {params.reference} {input.assembly} > {input.assembly}.paf 2>> {log}

            # Parse PAF and identify Ct contigs
            # PAF format: qname qlen qstart qend strand tname tlen tstart tend matches alnlen mapq ...
            # Keep contigs where: (matches/qlen)*100 >= min_coverage AND (matches/alnlen)*100 >= min_identity

            awk -v min_id={params.min_identity} -v min_cov={params.min_coverage} '
            BEGIN {{ OFS="\t" }}
            {{
                qname = $1
                qlen = $2
                matches = $10
                alnlen = $11

                # Calculate identity and coverage
                identity = (alnlen > 0) ? (matches / alnlen) * 100 : 0
                coverage = (qlen > 0) ? (alnlen / qlen) * 100 : 0

                # Track best alignment per contig
                if (!(qname in best_cov) || coverage > best_cov[qname]) {{
                    best_cov[qname] = coverage
                    best_id[qname] = identity
                    best_len[qname] = qlen
                }}
            }}
            END {{
                for (contig in best_cov) {{
                    if (best_id[contig] >= min_id && best_cov[contig] >= min_cov) {{
                        print contig, best_len[contig], best_id[contig], best_cov[contig], "KEEP"
                    }} else {{
                        print contig, best_len[contig], best_id[contig], best_cov[contig], "REMOVE"
                    }}
                }}
            }}' {input.assembly}.paf > {input.assembly}.filter_decisions.tsv

            # Also check for contigs with NO alignment (complete contamination)
            cut -f1 {input.assembly}.fai | while read contig; do
                if ! grep -q "^$contig" {input.assembly}.filter_decisions.tsv; then
                    len=$(grep "^$contig" {input.assembly}.fai | cut -f2)
                    echo -e "$contig\t$len\t0\t0\tREMOVE_NO_ALIGN"
                fi
            done >> {input.assembly}.filter_decisions.tsv

            # Extract contigs to keep (|| true prevents exit code 1 when no matches)
            grep "KEEP" {input.assembly}.filter_decisions.tsv | cut -f1 > {input.assembly}.keep_contigs.txt || true

            # Create filtered assembly
            if [ -s {input.assembly}.keep_contigs.txt ]; then
                seqtk subseq {input.assembly} {input.assembly}.keep_contigs.txt > {output.filtered}
            else
                # No contigs passed - keep original (shouldn't happen for real Ct samples)
                cp {input.assembly} {output.filtered}
                echo "WARNING: No contigs passed filtering - keeping original" >> {log}
            fi

            # Generate report
            total_contigs=$(wc -l < {input.assembly}.fai)
            kept_contigs=$(grep -c "KEEP" {input.assembly}.filter_decisions.tsv) || kept_contigs=0
            removed_contigs=$((total_contigs - kept_contigs))

            original_size=$(awk '{{sum+=$2}} END {{print sum}}' {input.assembly}.fai)
            filtered_size=$(seqtk comp {output.filtered} | awk '{{sum+=$2}} END {{print sum}}') || filtered_size=0
            removed_bp=$((original_size - filtered_size))

            cat > {output.report} << EOF
ASSEMBLY FILTERING REPORT
=========================
Sample: {params.sample}
Reference: {params.reference}
Min identity: {params.min_identity}%
Min coverage: {params.min_coverage}%

RESULTS
-------
Total contigs: $total_contigs
Kept contigs: $kept_contigs
Removed contigs: $removed_contigs

Original size: $original_size bp
Filtered size: $filtered_size bp
Removed: $removed_bp bp

CONTIG DETAILS
--------------
EOF

            echo -e "Contig\tLength\tIdentity\tCoverage\tStatus" >> {output.report}
            cat {input.assembly}.filter_decisions.tsv >> {output.report}

            # Cleanup temp files
            rm -f {input.assembly}.paf {input.assembly}.filter_decisions.tsv {input.assembly}.keep_contigs.txt {input.assembly}.fai

            touch {output.status}
        """


    rule apply_filtered_assembly:
        """
        Replace the original assembly with the filtered version.
        Keeps a backup of the unfiltered assembly.
        """
        input:
            filtered = OUTDIR / "{sample}" / "best" / "assembly.filtered.fasta",
            original = OUTDIR / "{sample}" / "best" / "assembly.fasta",
            report = OUTDIR / "{sample}" / "best" / "filter_report.txt"
        output:
            backup = OUTDIR / "{sample}" / "best" / "assembly.unfiltered.fasta",
            applied = OUTDIR / "status" / "filter_applied.{sample}.txt"
        shell: """
            # Backup original
            cp {input.original} {output.backup}

            # Replace with filtered
            cp {input.filtered} {input.original}

            touch {output.applied}
        """


    rule filter_summary:
        """Summary of filtering across all samples."""
        input:
            reports = expand(OUTDIR / "{sample}" / "best" / "filter_report.txt", sample=SAMPLES),
            applied = expand(OUTDIR / "status" / "filter_applied.{sample}.txt", sample=SAMPLES)
        output:
            summary = OUTDIR / "reports" / "assembly_filter_summary.tsv"
        run:
            import re

            rows = []
            for report_file in input.reports:
                sample = report_file.split("/")[-3]

                with open(report_file) as f:
                    content = f.read()

                # Parse values
                original = re.search(r"Original size: (\d+)", content)
                filtered = re.search(r"Filtered size: (\d+)", content)
                removed = re.search(r"Removed: (\d+)", content)
                kept = re.search(r"Kept contigs: (\d+)", content)
                total = re.search(r"Total contigs: (\d+)", content)

                rows.append({
                    "sample": sample,
                    "original_bp": original.group(1) if original else "NA",
                    "filtered_bp": filtered.group(1) if filtered else "NA",
                    "removed_bp": removed.group(1) if removed else "NA",
                    "total_contigs": total.group(1) if total else "NA",
                    "kept_contigs": kept.group(1) if kept else "NA"
                })

            with open(output.summary, 'w') as f:
                f.write("sample\toriginal_bp\tfiltered_bp\tremoved_bp\ttotal_contigs\tkept_contigs\n")
                for row in sorted(rows, key=lambda x: x["sample"]):
                    f.write(f"{row['sample']}\t{row['original_bp']}\t{row['filtered_bp']}\t{row['removed_bp']}\t{row['total_contigs']}\t{row['kept_contigs']}\n")

else:
    # Filtering disabled - create dummy status files
    rule skip_filter:
        input:
            assembly = OUTDIR / "{sample}" / "best" / "assembly.fasta"
        output:
            status = OUTDIR / "status" / "filter_applied.{sample}.txt"
        shell:
            "touch {output.status}"
