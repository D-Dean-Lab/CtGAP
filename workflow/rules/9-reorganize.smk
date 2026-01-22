# workflow/rules/9-reorganize.smk
# Final step: Create organized view of outputs using symlinks
# Original files stay in place (Snakemake-compatible), symlinks provide clean access

import os
from pathlib import Path

rule reorganize_outputs:
    """
    Create organized output structure using symlinks.

    This creates a clean user-facing view while keeping original files
    in place for Snakemake DAG compatibility.

    Structure created:
    output/
    ├── 00_README.txt
    ├── 01_per_sample/    → symlinks to sample dirs
    ├── 02_reports/       → symlinks to reports
    ├── 03_results/       → symlinks to collated results
    ├── 04_phylogeny/     → symlinks to tree files
    └── [original dirs remain for Snakemake]
    """
    input:
        # Wait for all main outputs to complete
        reports = expand(OUTDIR / "{sample}" / "reports" / "{sample}_ctgap_report.pdf", sample=SAMPLES),
        summary_csv = OUTDIR / "reports" / "ctgap_all_samples_summary.csv",
        batch_summary = OUTDIR / "reports" / "ctgap_batch_summary.txt",
        tree = OUTDIR / f"{ORG}.tree",
    output:
        status = OUTDIR / "status" / "reorganize_complete.txt",
        readme = OUTDIR / "00_README.txt"
    params:
        samples = SAMPLES,
        outdir = OUTDIR,
        org = ORG
    run:
        from pathlib import Path
        import os

        outdir = Path(params.outdir)

        # Define organized structure directories
        dirs = {
            "per_sample": outdir / "01_per_sample",
            "reports": outdir / "02_reports",
            "results": outdir / "03_results",
            "phylogeny": outdir / "04_phylogeny",
            "annotation": outdir / "05_annotation",
        }

        # Create directories
        for d in dirs.values():
            d.mkdir(parents=True, exist_ok=True)

        # Subdirectories
        (dirs["reports"] / "per_sample").mkdir(exist_ok=True)
        (dirs["reports"] / "summaries").mkdir(exist_ok=True)
        (dirs["reports"] / "qc").mkdir(exist_ok=True)
        (dirs["results"] / "typing").mkdir(exist_ok=True)
        (dirs["results"] / "coverage").mkdir(exist_ok=True)
        (dirs["results"] / "plasmid").mkdir(exist_ok=True)

        def make_symlink(src, dst):
            """Create relative symlink from dst pointing to src."""
            src = Path(src)
            dst = Path(dst)
            if not src.exists():
                return False
            dst.parent.mkdir(parents=True, exist_ok=True)
            if dst.exists() or dst.is_symlink():
                dst.unlink()
            try:
                rel_path = os.path.relpath(src, dst.parent)
                dst.symlink_to(rel_path)
                return True
            except Exception as e:
                print(f"Warning: Could not create symlink {dst} -> {src}: {e}")
                return False

        # =================================================================
        # 01_per_sample: Link sample directories
        # =================================================================
        for sample in params.samples:
            make_symlink(outdir / sample, dirs["per_sample"] / sample)

        # =================================================================
        # 02_reports: Organize report links
        # =================================================================

        # Per-sample comprehensive reports
        for sample in params.samples:
            make_symlink(
                outdir / sample / "reports" / f"{sample}_ctgap_report.pdf",
                dirs["reports"] / "per_sample" / f"{sample}_ctgap_report.pdf"
            )
            make_symlink(
                outdir / sample / "reports" / f"{sample}_plasmid_report.pdf",
                dirs["reports"] / "per_sample" / f"{sample}_plasmid_report.pdf"
            )

        # Summary reports
        summary_links = [
            ("reports/ctgap_all_samples_summary.csv", "summaries/ctgap_all_samples.csv"),
            ("reports/ctgap_batch_summary.txt", "summaries/batch_summary.txt"),
            ("reports/assembly_methods_summary.txt", "summaries/assembly_methods.txt"),
            ("reports/scaffold_summary.txt", "summaries/scaffold_summary.txt"),
            ("reports/ref_scaffold_summary.txt", "summaries/ref_scaffold_summary.txt"),
            ("reports/comprehensive_plasmid_analysis.txt", "summaries/plasmid_analysis.txt"),
            ("reports/ct_typing_all_samples.csv", "summaries/ct_typing_summary.csv"),
            ("reports/denovo.alignment_summary.txt", "summaries/denovo_alignment.txt"),
            ("reports/ref-denovo.alignment_summary.txt", "summaries/ref-denovo_alignment.txt"),
        ]
        for src_rel, dst_rel in summary_links:
            make_symlink(outdir / src_rel, dirs["reports"] / dst_rel)

        # Assembly selection
        make_symlink(
            outdir / "reports" / "assembly_selection",
            dirs["reports"] / "summaries" / "assembly_selection"
        )

        # QC reports
        make_symlink(outdir / "multiqc", dirs["reports"] / "qc" / "multiqc")

        # =================================================================
        # 03_results: Collated analysis results
        # =================================================================

        # Typing
        typing_links = [
            ("best.ompA_genovar.blast.tsv", "typing/all_samples.ompa_blast.tsv"),
            ("best.secondary.blast.tsv", "typing/all_samples.secondary_blast.tsv"),
            ("best.mlst.generic.results.tsv", "typing/all_samples.mlst_chlamydiales.tsv"),
            ("best.mlst.ct.results.tsv", "typing/all_samples.mlst_ctrachomatis.tsv"),
        ]
        for src_name, dst_rel in typing_links:
            make_symlink(outdir / src_name, dirs["results"] / dst_rel)

        # Coverage
        make_symlink(outdir / "denovo.coverage.tsv", dirs["results"] / "coverage" / "denovo.coverage.tsv")
        make_symlink(outdir / "ref-denovo.coverage.tsv", dirs["results"] / "coverage" / "ref-denovo.coverage.tsv")

        # Dedicated plasmid assembly (reference-guided + Unicycler)
        plasmid_links = [
            # Collated results
            ("reports/plasmid_summary.tsv", "plasmid/plasmid_summary.tsv"),
            ("plasmid.blast.all_samples.tsv", "plasmid/plasmid_blast.tsv"),
            ("plasmid.mlst.all_samples.tsv", "plasmid/plasmid_mlst.tsv"),
            ("reports/plasmid_analysis_report.txt", "plasmid/analysis_report.txt"),
        ]
        for src_name, dst_rel in plasmid_links:
            make_symlink(outdir / src_name, dirs["results"] / dst_rel)

        # Per-sample plasmid assemblies
        (dirs["results"] / "plasmid" / "per_sample").mkdir(exist_ok=True)
        for sample in params.samples:
            make_symlink(
                outdir / sample / "plasmid" / f"{sample}.plasmid.fasta",
                dirs["results"] / "plasmid" / "per_sample" / f"{sample}.plasmid.fasta"
            )

        # CT typing
        make_symlink(outdir / "ct_typing", dirs["results"] / "ct_typing")

        # =================================================================
        # 04_phylogeny: Tree files
        # =================================================================
        make_symlink(outdir / "tree", dirs["phylogeny"] / "tree_data")
        make_symlink(outdir / f"{params.org}.tree", dirs["phylogeny"] / f"{params.org}.tree")

        # =================================================================
        # 05_annotation: Bakta genome annotation (if enabled)
        # =================================================================
        # Link annotation summary reports
        make_symlink(outdir / "reports" / "annotation_summary.tsv", dirs["annotation"] / "annotation_summary.tsv")
        make_symlink(outdir / "reports" / "annotation_report.txt", dirs["annotation"] / "annotation_report.txt")
        make_symlink(outdir / "annotation" / "coverage_check.tsv", dirs["annotation"] / "coverage_check.tsv")

        # Per-sample annotation outputs (GFF3, GenBank, proteins)
        (dirs["annotation"] / "per_sample").mkdir(exist_ok=True)
        for sample in params.samples:
            sample_annotation_dir = outdir / sample / "annotation"
            if sample_annotation_dir.exists():
                make_symlink(sample_annotation_dir, dirs["annotation"] / "per_sample" / sample)

        # =================================================================
        # Mixed strain detection results (in reports)
        # =================================================================
        make_symlink(outdir / "reports" / "mixed_detection_summary.tsv", dirs["reports"] / "summaries" / "mixed_detection_summary.tsv")
        make_symlink(outdir / "reports" / "mixed_detection_report.txt", dirs["reports"] / "summaries" / "mixed_detection_report.txt")

        # =================================================================
        # Create README
        # =================================================================
        readme_content = f"""# CtGAP Output Directory
Generated by CtGAP Pipeline

## Organized Structure (Quick Access)

These directories contain symlinks to the actual data for easy navigation:

01_per_sample/     → Sample-specific outputs
   └── {{sample}}/
       ├── best/        - Selected best assembly + typing results
       ├── denovo/      - De novo assembly outputs
       ├── ref-denovo/  - Reference-guided assembly
       ├── ct_typing/   - CT typing results
       ├── qc/          - Quality control (fastp)
       ├── scrub/       - Host depletion (kraken)
       └── reports/     - Per-sample PDF reports

02_reports/        → All reports in one place
   ├── per_sample/ - Individual sample PDFs
   ├── summaries/  - Cross-sample summary tables
   └── qc/         - MultiQC report

03_results/        → Collated analysis results
   ├── typing/     - ompA BLAST, MLST results (all samples)
   ├── coverage/   - Coverage statistics
   └── plasmid/    - Plasmid BLAST, MLST results

04_phylogeny/      → Phylogenetic analysis
   ├── {params.org}.tree    - Final tree (Newick format)
   └── tree_data/           - SKA alignment and inputs

## Quick Access Commands

# View all sample reports
ls 02_reports/per_sample/

# Open summary table
cat 02_reports/summaries/ctgap_all_samples.csv

# View genovar assignments
cat 03_results/typing/all_samples.ompa_blast.tsv

# View MLST results
cat 03_results/typing/all_samples.mlst_ctrachomatis.tsv

# Copy tree for visualization
cp 04_phylogeny/{params.org}.tree ~/my_tree.nwk

## Original Data Locations

The original output files remain in their standard locations for
Snakemake compatibility. The 01-04 directories are symlinks only.

## Samples Processed: {len(list(params.samples))}

{chr(10).join(f"  - {s}" for s in sorted(params.samples))}

## Pipeline Info

Mode: {shell("grep 'mode:' config/config.yaml | head -1", read=True).strip() if False else "See config.yaml"}
Date: {shell("date", read=True).strip() if False else "See log files"}
"""
        # Simplified readme without shell calls
        readme_content = f"""# CtGAP Output Directory
Generated by CtGAP Pipeline

## Organized Structure (Quick Access)

01_per_sample/     → Sample-specific outputs
02_reports/        → All reports in one place
03_results/        → Collated analysis results
04_phylogeny/      → Phylogenetic analysis
05_annotation/     → Genome annotations (Bakta)

## Key Files

- Summary table:    02_reports/summaries/ctgap_all_samples.csv
- Batch summary:    02_reports/summaries/batch_summary.txt
- Genovar results:  03_results/typing/all_samples.ompa_blast.tsv
- MLST results:     03_results/typing/all_samples.mlst_ctrachomatis.tsv
- Mixed detection:  02_reports/summaries/mixed_detection_summary.tsv
- Phylogeny:        04_phylogeny/{params.org}.tree
- Annotations:      05_annotation/annotation_summary.tsv

## Per-Sample Reports

PDF reports for each sample are in: 02_reports/per_sample/

## Genome Annotation

Bakta annotations are generated for samples with >=90% genome coverage.
- GFF3, GenBank, and protein FASTA files: 05_annotation/per_sample/
- Summary: 05_annotation/annotation_report.txt

## Mixed Strain Detection

Samples are analyzed for potential co-infections via variant heterozygosity.
- Samples with >50 heterozygous sites are flagged as potential mixed infections
- Check: 02_reports/summaries/mixed_detection_report.txt
- Per-sample reports show prominent warnings for flagged samples

## Samples Processed: {len(list(params.samples))}

{chr(10).join(f"  - {s}" for s in sorted(params.samples))}

## Notes

- Directories 01-05 contain symlinks to actual data
- Original files remain in place for Snakemake compatibility
- All data is accessible through either path structure
"""
        with open(output.readme, "w") as f:
            f.write(readme_content)

        # Touch status
        Path(output.status).touch()
