# workflow/rules/8-reports.smk
# Comprehensive reporting - combines chromosome and plasmid analysis

# -----------------------------------------------------------------------------
# 1. COMPREHENSIVE PER-SAMPLE PDF REPORT
# -----------------------------------------------------------------------------

rule create_comprehensive_report:
    """
    Generate a comprehensive PDF report for each sample combining:
    - Assembly selection and metrics
    - ompA genotype typing
    - MLST results
    - Ct strain typing results
    - Plasmid analysis
    - QC metrics
    """
    input:
        # Assembly
        best_assembly = OUTDIR / "{sample}" / "best" / "assembly.fasta",
        assembly_method = OUTDIR / "reports" / "assembly_selection" / "{sample}.method.txt",
        assembly_report = OUTDIR / "reports" / "assembly_selection" / "{sample}.txt",
        # QUAST - get from whichever mode produced the best
        denovo_quast = OUTDIR / "{sample}" / "denovo" / "assembly_statistics" / "report.tsv",
        # Typing
        blast_ompa = OUTDIR / "{sample}" / "best" / "blast" / "blast.ompa.tab",
        blast_secondary = OUTDIR / "{sample}" / "best" / "blast" / "blast.secondary.tab",
        mlst_chlamydiales = OUTDIR / "{sample}" / "best" / "mlst" / "{sample}.genome.chlamydiales.mlst.txt",
        mlst_ct = OUTDIR / "{sample}" / "best" / "mlst" / "{sample}.genome.ctrachomatis.mlst.txt",
        # Ct strain typing (if enabled)
        ct_typing = OUTDIR / "{sample}" / "ct_typing" / "{sample}.ct_assignment.csv" if config.get("ct_typing", {}).get("enabled", False) else [],
        # Plasmid (dedicated reference-guided assembly)
        plasmid_analysis = OUTDIR / "{sample}" / "plasmid" / "{sample}.plasmid_analysis.txt",
        plasmid_blast = OUTDIR / "{sample}" / "plasmid" / "{sample}.plasmid_blast.tsv",
        plasmid_mlst = OUTDIR / "{sample}" / "plasmid" / "{sample}.plasmid_mlst.txt",
        # QC
        fastp_json = OUTDIR / "{sample}" / "qc" / "{sample}_scrub.json",
        kraken_report = OUTDIR / "{sample}" / "scrub" / "{sample}.kraken2.report",
        # Coverage
        coverage = OUTDIR / "{sample}" / "denovo" / "bowtie_ref24" / "{sample}.coverage.ref24.tsv",
        # Annotation completion (wait for annotation to finish so we can include results)
        annotation_complete = OUTDIR / "status" / "annotation.complete.txt" if config.get("annotation", {}).get("enabled", True) else [],
        # Mixed detection completion (wait for mixed strain analysis)
        mixed_complete = OUTDIR / "status" / "mixed_detection.complete.txt" if config.get("mixed_detection", {}).get("enabled", True) else [],
    output:
        pdf = OUTDIR / "{sample}" / "reports" / "{sample}_ctgap_report.pdf",
        status = OUTDIR / "status" / "comprehensive_report.{sample}.txt"
    params:
        sample = lambda w: w.sample,
        ct_typing_enabled = config.get("ct_typing", {}).get("enabled", False),
        annotation_enabled = config.get("annotation", {}).get("enabled", True),
        mixed_detection_enabled = config.get("mixed_detection", {}).get("enabled", True),
        outdir = str(OUTDIR)
    conda: "../envs/pdf.yaml"
    log: OUTDIR / "{sample}" / "log" / "comprehensive_report.{sample}.log"
    script: "../scripts/create_comprehensive_report.py"


# -----------------------------------------------------------------------------
# 2. ALL-SAMPLES SUMMARY CSV
# -----------------------------------------------------------------------------

rule create_summary_csv:
    """
    Generate a master CSV summarizing all samples with key metrics.
    One row per sample for easy comparison and downstream analysis.
    """
    input:
        # Wait for all individual reports
        reports = expand(OUTDIR / "{sample}" / "reports" / "{sample}_ctgap_report.pdf", sample=SAMPLES),
        # Assembly methods
        methods = expand(OUTDIR / "reports" / "assembly_selection" / "{sample}.method.txt", sample=SAMPLES),
        # Assembly reports (contains scores)
        assembly_reports = expand(OUTDIR / "reports" / "assembly_selection" / "{sample}.txt", sample=SAMPLES),
        # BLAST results
        blast_collated = OUTDIR / "best.ompA_genovar.blast.tsv",
        # MLST results
        mlst_chlamydiales = OUTDIR / "best.mlst.generic.results.tsv",
        mlst_ct = OUTDIR / "best.mlst.ct.results.tsv",
        # Ct strain typing (if enabled)
        ct_typing_summary = OUTDIR / "reports" / "ct_typing_all_samples.csv" if config.get("ct_typing", {}).get("enabled", False) else [],
        # QUAST reports (for assembly metrics)
        quast_reports = expand(OUTDIR / "{sample}" / "denovo" / "assembly_statistics" / "report.tsv", sample=SAMPLES),
        # Plasmid analysis (dedicated reference-guided assembly)
        plasmid_analyses = expand(OUTDIR / "{sample}" / "plasmid" / "{sample}.plasmid_analysis.txt", sample=SAMPLES),
    output:
        csv = OUTDIR / "reports" / "ctgap_all_samples_summary.csv",
        status = OUTDIR / "status" / "summary_csv.txt"
    params:
        samples = SAMPLES,
        outdir = OUTDIR,
        ct_typing_enabled = config.get("ct_typing", {}).get("enabled", False)
    log: OUTDIR / "log" / "summary_csv.log"
    run:
        import csv
        import os
        import re
        from pathlib import Path

        rows = []

        for sample in params.samples:
            row = {"sample": sample}

            # 1. Assembly method
            method_file = Path(params.outdir) / "reports" / "assembly_selection" / f"{sample}.method.txt"
            if method_file.exists():
                row["assembly_method"] = method_file.read_text().strip()
            else:
                row["assembly_method"] = "unknown"

            # 2. Assembly score (parse from report)
            report_file = Path(params.outdir) / "reports" / "assembly_selection" / f"{sample}.txt"
            row["assembly_score"] = ""
            if report_file.exists():
                content = report_file.read_text()
                # Look for score in report
                match = re.search(r"Score:\s*([\d.]+)", content)
                if match:
                    row["assembly_score"] = match.group(1)

            # 3. QUAST metrics
            quast_file = Path(params.outdir) / sample / "denovo" / "assembly_statistics" / "report.tsv"
            row["total_length"] = ""
            row["num_contigs"] = ""
            row["n50"] = ""
            row["gc_percent"] = ""
            if quast_file.exists():
                with open(quast_file) as f:
                    for line in f:
                        parts = line.strip().split("\t")
                        if len(parts) >= 2:
                            key, val = parts[0], parts[1]
                            if key == "Total length":
                                row["total_length"] = val.replace(",", "")
                            elif key == "# contigs":
                                row["num_contigs"] = val
                            elif key == "N50":
                                row["n50"] = val.replace(",", "")
                            elif key == "GC (%)":
                                row["gc_percent"] = val

            # 4. ompA genotype (parse from BLAST)
            blast_file = Path(params.outdir) / sample / "best" / "blast" / "blast.ompa.tab"
            row["genotype"] = ""
            row["genotype_identity"] = ""
            if blast_file.exists() and blast_file.stat().st_size > 0:
                with open(blast_file) as f:
                    lines = f.readlines()
                    if lines:
                        # Get top hit by bitscore
                        best_line = sorted(lines, key=lambda x: float(x.split("\t")[11]) if len(x.split("\t")) > 11 else 0, reverse=True)[0]
                        fields = best_line.strip().split("\t")
                        if len(fields) >= 3:
                            # Extract genotype from subject ID (e.g., "E_AF063199.1" -> "E")
                            subject = fields[1]
                            genotype_match = re.match(r"([A-Za-z]+)_", subject)
                            if genotype_match:
                                row["genotype"] = genotype_match.group(1)
                            row["genotype_identity"] = fields[2]

            # 5. MLST (simplified - just ST)
            mlst_file = Path(params.outdir) / sample / "best" / "mlst" / f"{sample}.genome.ctrachomatis.mlst.txt"
            row["mlst_st"] = ""
            if mlst_file.exists() and mlst_file.stat().st_size > 0:
                # claMLST output is TSV with header row
                with open(mlst_file) as mf:
                    lines = [l.strip() for l in mf.readlines() if l.strip()]
                    # Skip header row, get data from second line
                    if len(lines) > 1:
                        data_line = lines[1]
                        fields = data_line.split("\t")
                        # claMLST format: Sample, ST, allele1, allele2, ...
                        if len(fields) >= 2 and fields[1] != "ST":
                            st_val = fields[1]
                            if st_val and st_val != "-" and st_val != "NA":
                                row["mlst_st"] = f"ST{st_val}" if not st_val.startswith("ST") else st_val

            # 6. Ct genome backbone/strain typing (if enabled)
            row["genome_backbone"] = ""
            row["backbone_confidence"] = ""
            if params.ct_typing_enabled:
                ct_file = Path(params.outdir) / sample / "ct_typing" / f"{sample}.ct_assignment.csv"
                if ct_file.exists() and ct_file.stat().st_size > 0:
                    # Parse Ct assignment CSV (comma-delimited)
                    import csv as csv_mod
                    with open(ct_file) as cf:
                        reader = csv_mod.DictReader(cf)
                        for ct_row in reader:
                            # Get the assigned backbone/strain from putative_reference column
                            if "putative_reference" in ct_row:
                                row["genome_backbone"] = ct_row.get("putative_reference", "")
                            elif "ref_clean" in ct_row:
                                row["genome_backbone"] = ct_row.get("ref_clean", "")
                            row["backbone_confidence"] = ct_row.get("confidence", "")
                            break  # Only first row

            # 7. Plasmid detection status
            plasmid_file = Path(params.outdir) / sample / "plasmid" / f"{sample}.plasmid_analysis.txt"
            row["plasmid_status"] = "Unknown"
            row["plasmid_length"] = ""
            if plasmid_file.exists():
                content = plasmid_file.read_text()
                if "CT PLASMID DETECTED" in content:
                    row["plasmid_status"] = "DETECTED"
                    # Extract length
                    for line in content.split("\n"):
                        if "Length:" in line:
                            match = re.search(r"Length:\s*(\d+)", line)
                            if match:
                                row["plasmid_length"] = match.group(1)
                            break
                elif "NO PLASMID DETECTED" in content:
                    row["plasmid_status"] = "NOT_DETECTED"
                elif "NO VALID CT PLASMID" in content:
                    row["plasmid_status"] = "NO_VALID"

            rows.append(row)

        # Write CSV
        fieldnames = [
            "sample", "genotype", "genotype_identity", "mlst_st",
            "genome_backbone", "backbone_confidence", "assembly_method", "assembly_score",
            "total_length", "num_contigs", "n50", "gc_percent",
            "plasmid_status", "plasmid_length"
        ]

        with open(output.csv, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction='ignore')
            writer.writeheader()
            writer.writerows(rows)

        # Touch status file
        Path(output.status).touch()


# -----------------------------------------------------------------------------
# 3. BATCH SUMMARY REPORT (TEXT)
# -----------------------------------------------------------------------------

rule create_batch_summary:
    """
    Generate a human-readable batch summary report.
    """
    input:
        csv = OUTDIR / "reports" / "ctgap_all_samples_summary.csv"
    output:
        report = OUTDIR / "reports" / "ctgap_batch_summary.txt"
    run:
        import csv
        from collections import Counter
        from datetime import datetime

        # Read CSV
        with open(input.csv) as f:
            reader = csv.DictReader(f)
            rows = list(reader)

        # Generate summary
        with open(output.report, "w") as out:
            out.write("=" * 70 + "\n")
            out.write("CtGAP BATCH ANALYSIS SUMMARY\n")
            out.write("=" * 70 + "\n")
            out.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            out.write(f"Total samples: {len(rows)}\n")
            out.write("\n")

            # Genotype distribution
            genotypes = Counter(r.get("genotype", "Unknown") for r in rows if r.get("genotype"))
            out.write("-" * 40 + "\n")
            out.write("GENOTYPE DISTRIBUTION\n")
            out.write("-" * 40 + "\n")
            for geno, count in sorted(genotypes.items()):
                out.write(f"  {geno}: {count} samples ({100*count/len(rows):.1f}%)\n")
            out.write("\n")

            # Assembly method distribution
            methods = Counter(r.get("assembly_method", "unknown") for r in rows)
            out.write("-" * 40 + "\n")
            out.write("ASSEMBLY METHOD SELECTION\n")
            out.write("-" * 40 + "\n")
            for method, count in sorted(methods.items()):
                out.write(f"  {method}: {count} samples ({100*count/len(rows):.1f}%)\n")
            out.write("\n")

            # Ct strain typing (if available)
            strains = Counter(r.get("genome_backbone", "") for r in rows if r.get("genome_backbone"))
            if strains:
                out.write("-" * 40 + "\n")
                out.write("Ct STRAIN DISTRIBUTION\n")
                out.write("-" * 40 + "\n")
                for strain, count in sorted(strains.items()):
                    out.write(f"  {strain}: {count} samples\n")
                out.write("\n")

            # Plasmid summary
            plasmid_status = Counter(r.get("plasmid_status", "Unknown") for r in rows)
            out.write("-" * 40 + "\n")
            out.write("PLASMID DETECTION\n")
            out.write("-" * 40 + "\n")
            for status, num in sorted(plasmid_status.items()):
                out.write(f"  {status}: {num} samples ({100*num/len(rows):.1f}%)\n")
            out.write("\n")

            # Assembly quality summary
            out.write("-" * 40 + "\n")
            out.write("ASSEMBLY QUALITY METRICS\n")
            out.write("-" * 40 + "\n")
            lengths = [int(r["total_length"]) for r in rows if r.get("total_length", "").isdigit()]
            n50s = [int(r["n50"]) for r in rows if r.get("n50", "").isdigit()]
            contigs = [int(r["num_contigs"]) for r in rows if r.get("num_contigs", "").isdigit()]

            if lengths:
                out.write(f"  Total length: {min(lengths):,} - {max(lengths):,} bp (median: {sorted(lengths)[len(lengths)//2]:,})\n")
            if n50s:
                out.write(f"  N50: {min(n50s):,} - {max(n50s):,} bp (median: {sorted(n50s)[len(n50s)//2]:,})\n")
            if contigs:
                out.write(f"  Contigs: {min(contigs)} - {max(contigs)} (median: {sorted(contigs)[len(contigs)//2]})\n")
            out.write("\n")

            out.write("=" * 70 + "\n")
            out.write("For detailed per-sample reports, see: reports/per_sample/\n")
            out.write("For full data, see: reports/ctgap_all_samples_summary.csv\n")
            out.write("=" * 70 + "\n")
