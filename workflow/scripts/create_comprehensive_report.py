#!/usr/bin/env python3
"""
Create comprehensive per-sample PDF report for CtGAP.
Combines chromosome and plasmid analysis into a single document.
"""

import os
import sys
import json
import re
from pathlib import Path
from datetime import datetime

from reportlab.lib.pagesizes import letter, A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_LEFT, TA_RIGHT
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle,
    PageBreak, HRFlowable, KeepTogether
)


def parse_quast_report(quast_file):
    """Parse QUAST TSV report."""
    stats = {}
    if not os.path.exists(quast_file):
        return stats

    with open(quast_file) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                key, val = parts[0], parts[1]
                if key == "Total length":
                    stats["total_length"] = int(val.replace(",", ""))
                elif key == "# contigs":
                    stats["num_contigs"] = int(val)
                elif key == "N50":
                    stats["n50"] = int(val.replace(",", ""))
                elif key == "GC (%)":
                    stats["gc_percent"] = float(val)
                elif key == "Largest contig":
                    stats["largest_contig"] = int(val.replace(",", ""))
    return stats


def count_gaps_in_main_contig(assembly_file):
    """Count N's (gaps) in the largest contig of the assembly."""
    if not os.path.exists(assembly_file):
        return 0

    # Parse FASTA and find largest contig
    largest_seq = ""
    current_seq = ""
    current_header = ""

    with open(assembly_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                # Save previous contig if larger
                if len(current_seq) > len(largest_seq):
                    largest_seq = current_seq
                current_header = line[1:]
                current_seq = ""
            else:
                current_seq += line

        # Check last contig
        if len(current_seq) > len(largest_seq):
            largest_seq = current_seq

    # Count N's in largest contig
    return largest_seq.upper().count('N')


def parse_blast_results(blast_file):
    """Parse BLAST tabular output."""
    results = []
    if not os.path.exists(blast_file) or os.path.getsize(blast_file) == 0:
        return results

    with open(blast_file) as f:
        for line in f:
            fields = line.strip().split("\t")
            if len(fields) >= 12:
                # Skip header lines (check if identity field is numeric)
                try:
                    identity = float(fields[2])
                except ValueError:
                    continue  # Skip header or malformed lines

                results.append({
                    "query": fields[0],
                    "subject": fields[1],
                    "identity": identity,
                    "length": int(fields[3]),
                    "evalue": float(fields[10]),
                    "bitscore": float(fields[11])
                })

    # Sort by bitscore descending
    results.sort(key=lambda x: x["bitscore"], reverse=True)
    return results


def parse_mlst_results(mlst_file):
    """Parse claMLST output (TSV format with header)."""
    if not os.path.exists(mlst_file) or os.path.getsize(mlst_file) == 0:
        return "Not available"

    with open(mlst_file) as f:
        lines = [l.strip() for l in f.readlines() if l.strip()]

    # claMLST output has header row, data in second row
    # Format: Sample, ST, allele1, allele2, ...
    if len(lines) > 1:
        data_line = lines[1]
        fields = data_line.split("\t")
        if len(fields) >= 2 and fields[1] != "ST":
            st_val = fields[1]
            if st_val and st_val != "-" and st_val != "NA":
                return f"ST{st_val}" if not st_val.startswith("ST") else st_val
            else:
                return "No ST match"

    return "No match found"


def parse_assembly_selection(report_file, method_file):
    """Parse assembly selection report."""
    result = {"method": "unknown", "score": "N/A"}

    if os.path.exists(method_file):
        result["method"] = Path(method_file).read_text().strip()

    if os.path.exists(report_file):
        content = Path(report_file).read_text()
        # Extract scores (non-greedy .*? to match first Score after each header)
        denovo_match = re.search(r"De Novo.*?Score:\s*([\d.]+)", content, re.DOTALL)
        ref_match = re.search(r"Reference-Guided.*?Score:\s*([\d.]+)", content, re.DOTALL)

        if denovo_match:
            result["denovo_score"] = denovo_match.group(1)
        if ref_match:
            result["ref_score"] = ref_match.group(1)

        # Get winning score
        score_match = re.search(r"SELECTED:.*Score:\s*([\d.]+)", content)
        if score_match:
            result["score"] = score_match.group(1)

    return result


def parse_ct_typing(ct_file):
    """Parse CT typing assignment CSV."""
    if not os.path.exists(ct_file) or os.path.getsize(ct_file) == 0:
        return None

    result = {}

    try:
        import csv
        with open(ct_file) as f:
            # CT typing CSV is comma-delimited (default pandas to_csv)
            reader = csv.DictReader(f)
            for row in reader:
                # Get type from putative_reference column (primary)
                if "putative_reference" in row and row["putative_reference"]:
                    result["type"] = row["putative_reference"]
                elif "ref_clean" in row and row["ref_clean"]:
                    result["type"] = row["ref_clean"]

                # Get confidence
                if "confidence" in row and row["confidence"]:
                    result["confidence"] = row["confidence"]

                # Get ANI value
                if "ANI" in row and row["ANI"]:
                    try:
                        result["fastani_best"] = f"{float(row['ANI']):.2f}%"
                    except ValueError:
                        result["fastani_best"] = row["ANI"]

                # Get Mash distance
                if "mash_dist" in row and row["mash_dist"]:
                    result["mash_best"] = row["mash_dist"]

                break  # Only first row
    except Exception as e:
        print(f"Error parsing CT typing file: {e}", file=sys.stderr)

    return result if result else None


def parse_plasmid_analysis(analysis_file):
    """Parse dedicated plasmid analysis report."""
    result = {
        "status": "Unknown",
        "plasmid_count": 0,
        "length": "",
        "circular": "",
        "deviation": ""
    }

    if not os.path.exists(analysis_file):
        return result

    content = Path(analysis_file).read_text()

    # Determine status
    if "CT PLASMID DETECTED" in content:
        result["status"] = "DETECTED"
        result["plasmid_count"] = 1
    elif "NO PLASMID DETECTED" in content:
        result["status"] = "NOT_DETECTED"
    elif "NO VALID CT PLASMID" in content:
        result["status"] = "NO_VALID"

    # Extract details
    for line in content.split("\n"):
        if "Length:" in line and "bp" in line:
            match = re.search(r"Length:\s*(\d+)", line)
            if match:
                result["length"] = match.group(1)
        if "Circular:" in line:
            result["circular"] = line.split(":")[1].strip()
        if "Size deviation:" in line:
            result["deviation"] = line.split(":")[1].strip()

    return result


def parse_fastp_json(json_file):
    """Parse fastp JSON output for QC metrics."""
    result = {}

    if not os.path.exists(json_file):
        return result

    try:
        with open(json_file) as f:
            data = json.load(f)

        if "summary" in data:
            before = data["summary"].get("before_filtering", {})
            after = data["summary"].get("after_filtering", {})

            result["raw_reads"] = before.get("total_reads", 0)
            result["filtered_reads"] = after.get("total_reads", 0)
            result["raw_bases"] = before.get("total_bases", 0)
            result["filtered_bases"] = after.get("total_bases", 0)
            result["q30_rate"] = after.get("q30_rate", 0)
            result["gc_content"] = after.get("gc_content", 0)
    except Exception as e:
        print(f"Error parsing fastp JSON: {e}", file=sys.stderr)

    return result


def parse_kraken_report(report_file):
    """Parse Kraken2 report for top classifications."""
    results = []

    if not os.path.exists(report_file):
        return results

    with open(report_file) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 6:
                pct = float(parts[0])
                name = parts[5].strip()
                if pct >= 1.0 and name:  # Only significant hits
                    results.append({"name": name, "percent": pct})

    return sorted(results, key=lambda x: x["percent"], reverse=True)[:10]


def parse_coverage(coverage_file):
    """Parse samtools coverage output."""
    result = {}

    if not os.path.exists(coverage_file):
        return result

    with open(coverage_file) as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 7:
                result["mean_depth"] = float(parts[6]) if parts[6] != "nan" else 0
                result["coverage_pct"] = float(parts[5]) if parts[5] != "nan" else 0
                break

    return result


def parse_bakta_annotation(annotation_txt):
    """Parse Bakta annotation summary (.txt file)."""
    result = {}

    if not os.path.exists(annotation_txt):
        return None

    try:
        content = Path(annotation_txt).read_text()

        # Parse Bakta summary statistics
        patterns = {
            "sequences": r"Sequences:\s*(\d+)",
            "size": r"Size:\s*([\d,]+)",
            "gc": r"GC:\s*([\d.]+)",
            "cds": r"CDSs:\s*(\d+)",
            "trnas": r"tRNAs:\s*(\d+)",
            "rrnas": r"rRNAs:\s*(\d+)",
            "ncrnas": r"ncRNAs:\s*(\d+)",
            "tmrnas": r"tmRNAs:\s*(\d+)",
            "crispr": r"CRISPRs:\s*(\d+)",
            "hypotheticals": r"hypotheticals:\s*(\d+)",
            "pseudogenes": r"pseudogenes:\s*(\d+)",
        }

        for key, pattern in patterns.items():
            match = re.search(pattern, content, re.IGNORECASE)
            if match:
                val = match.group(1).replace(",", "")
                result[key] = val

        return result if result else None
    except Exception as e:
        print(f"Error parsing Bakta annotation: {e}", file=sys.stderr)
        return None


def parse_mixed_detection(summary_file):
    """Parse mixed strain detection summary."""
    result = {}

    if not os.path.exists(summary_file):
        return None

    try:
        content = Path(summary_file).read_text()

        # Extract heterozygous site count
        match = re.search(r"Heterozygous-like sites:\s*(\d+)", content)
        if match:
            result["het_sites"] = int(match.group(1))

        # Extract mean allele frequency
        match = re.search(r"Mean allele frequency:\s*([\d.]+|N/A)", content)
        if match:
            result["mean_af"] = match.group(1)

        # Extract main contig info (new format: "Main contig analyzed: contig_name (size bp)")
        match = re.search(r"Main contig analyzed:\s*(\S+)\s*\((\d+)\s*bp\)", content)
        if match:
            result["main_contig_name"] = match.group(1)
            result["main_contig_size"] = int(match.group(2))

        # Extract position range
        match = re.search(r"Position range:\s*(\d+)", content)
        if match:
            result["position_range"] = int(match.group(1))

        # Extract distribution (GENOME-WIDE or LOCALIZED)
        if "GENOME-WIDE" in content:
            result["distribution"] = "GENOME-WIDE"
        elif "LOCALIZED" in content:
            result["distribution"] = "LOCALIZED"

        # Extract status
        if "STATUS: MIXED_STRAIN_WARNING" in content:
            result["status"] = "MIXED_WARNING"
        elif "STATUS: LOW_HETEROZYGOSITY" in content:
            result["status"] = "LOW_HET"
        elif "STATUS: SINGLE_STRAIN" in content:
            result["status"] = "SINGLE"
        else:
            result["status"] = "Unknown"

        return result if result else None
    except Exception as e:
        print(f"Error parsing mixed detection: {e}", file=sys.stderr)
        return None


def create_pdf_report(sample, output_pdf, data):
    """Create the PDF report."""

    doc = SimpleDocTemplate(
        output_pdf,
        pagesize=letter,
        rightMargin=0.75*inch,
        leftMargin=0.75*inch,
        topMargin=0.75*inch,
        bottomMargin=0.75*inch
    )

    styles = getSampleStyleSheet()

    # Custom styles
    title_style = ParagraphStyle(
        'Title',
        parent=styles['Heading1'],
        fontSize=18,
        alignment=TA_CENTER,
        spaceAfter=12
    )

    section_style = ParagraphStyle(
        'Section',
        parent=styles['Heading2'],
        fontSize=14,
        textColor=colors.HexColor("#2C3E50"),
        spaceBefore=12,
        spaceAfter=6
    )

    subsection_style = ParagraphStyle(
        'Subsection',
        parent=styles['Heading3'],
        fontSize=11,
        textColor=colors.HexColor("#34495E"),
        spaceBefore=8,
        spaceAfter=4
    )

    normal_style = styles['Normal']

    key_value_style = ParagraphStyle(
        'KeyValue',
        parent=normal_style,
        fontSize=10,
        leftIndent=20
    )

    elements = []

    # =========================================================================
    # TITLE PAGE / HEADER
    # =========================================================================
    elements.append(Paragraph("CtGAP Analysis Report", title_style))
    elements.append(Spacer(1, 0.1*inch))

    # Header info table
    header_data = [
        ["Sample ID:", sample],
        ["Report Date:", datetime.now().strftime("%Y-%m-%d %H:%M")],
        ["Pipeline:", "CtGAP v2.0"],
    ]
    header_table = Table(header_data, colWidths=[1.5*inch, 4*inch])
    header_table.setStyle(TableStyle([
        ('FONTNAME', (0, 0), (0, -1), 'Helvetica-Bold'),
        ('FONTSIZE', (0, 0), (-1, -1), 10),
        ('BOTTOMPADDING', (0, 0), (-1, -1), 4),
        ('TOPPADDING', (0, 0), (-1, -1), 4),
    ]))
    elements.append(header_table)
    elements.append(Spacer(1, 0.2*inch))

    # Divider
    elements.append(HRFlowable(width="100%", thickness=2, color=colors.HexColor("#3498DB")))
    elements.append(Spacer(1, 0.2*inch))

    # =========================================================================
    # KEY FINDINGS SUMMARY BOX
    # =========================================================================
    elements.append(Paragraph("Key Findings", section_style))

    assembly = data.get("assembly_selection", {})
    blast = data.get("blast_ompa", [])
    ct = data.get("ct_typing")
    plasmid = data.get("plasmid_stats", {})
    quast = data.get("quast", {})

    # Extract genovar from BLAST
    genovar = "Unknown"
    genovar_identity = "N/A"
    if blast:
        subject = blast[0].get("subject", "")
        match = re.match(r"([A-Za-z]+)_", subject)
        if match:
            genovar = match.group(1)
        genovar_identity = f"{blast[0].get('identity', 0):.1f}%"

    # Plasmid status summary
    plasmid_summary = plasmid.get("status", "Unknown")
    if plasmid_summary == "DETECTED" and plasmid.get("length"):
        plasmid_summary = f"DETECTED ({plasmid.get('length')} bp)"

    findings_data = [
        ["ompA Genotype:", genovar, "Identity:", genovar_identity],
        ["Ct Strain:", ct.get("type", "N/A") if ct else "Not run", "Assembly Method:", assembly.get("method", "unknown")],
        ["Plasmid:", plasmid_summary, "Contigs:", str(quast.get("num_contigs", "N/A"))],
    ]

    findings_table = Table(findings_data, colWidths=[1.3*inch, 1.7*inch, 1.5*inch, 1.7*inch])
    findings_table.setStyle(TableStyle([
        ('FONTNAME', (0, 0), (0, -1), 'Helvetica-Bold'),
        ('FONTNAME', (2, 0), (2, -1), 'Helvetica-Bold'),
        ('FONTSIZE', (0, 0), (-1, -1), 10),
        ('BACKGROUND', (0, 0), (-1, -1), colors.HexColor("#ECF0F1")),
        ('BOX', (0, 0), (-1, -1), 1, colors.HexColor("#BDC3C7")),
        ('INNERGRID', (0, 0), (-1, -1), 0.5, colors.HexColor("#BDC3C7")),
        ('BOTTOMPADDING', (0, 0), (-1, -1), 6),
        ('TOPPADDING', (0, 0), (-1, -1), 6),
        ('LEFTPADDING', (0, 0), (-1, -1), 8),
    ]))
    elements.append(findings_table)
    elements.append(Spacer(1, 0.3*inch))

    # =========================================================================
    # SECTION 1: ASSEMBLY ANALYSIS
    # =========================================================================
    elements.append(Paragraph("1. Assembly Analysis", section_style))

    # Assembly selection
    elements.append(Paragraph("Assembly Selection", subsection_style))

    selection_text = f"""
    <b>Selected Method:</b> {assembly.get('method', 'unknown').replace('-', ' ').title()}
    """
    elements.append(Paragraph(selection_text, key_value_style))
    elements.append(Spacer(1, 0.1*inch))

    # QUAST metrics
    elements.append(Paragraph("Assembly Quality Metrics (QUAST)", subsection_style))

    if quast:
        # Get gap count from data (calculated from assembly FASTA)
        gap_count = quast.get('main_contig_gaps', 0)

        quast_data = [
            ["Metric", "Value"],
            ["Total Length", f"{quast.get('total_length', 0):,} bp"],
            ["Number of Contigs", str(quast.get('num_contigs', 'N/A'))],
            ["N50", f"{quast.get('n50', 0):,} bp"],
            ["Gaps in Main Contig", f"{gap_count:,} bp"],
            ["GC Content", f"{quast.get('gc_percent', 0):.2f}%"],
        ]
        quast_table = Table(quast_data, colWidths=[2.5*inch, 2.5*inch])
        quast_table.setStyle(TableStyle([
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor("#3498DB")),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.white),
            ('FONTSIZE', (0, 0), (-1, -1), 9),
            ('BOX', (0, 0), (-1, -1), 1, colors.HexColor("#BDC3C7")),
            ('INNERGRID', (0, 0), (-1, -1), 0.5, colors.HexColor("#BDC3C7")),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 4),
            ('TOPPADDING', (0, 0), (-1, -1), 4),
            ('ALIGN', (1, 1), (1, -1), 'RIGHT'),
        ]))
        elements.append(quast_table)
    else:
        elements.append(Paragraph("QUAST metrics not available", key_value_style))

    elements.append(Spacer(1, 0.2*inch))

    # =========================================================================
    # SECTION 2: TYPING RESULTS
    # =========================================================================
    elements.append(Paragraph("2. Typing Results", section_style))

    # ompA Genotype
    elements.append(Paragraph("ompA Genotype (BLAST)", subsection_style))

    if blast:
        top_hit = blast[0]
        blast_text = f"""
        <b>Best Match:</b> {top_hit.get('subject', 'N/A')}<br/>
        <b>Identity:</b> {top_hit.get('identity', 0):.2f}%<br/>
        <b>Alignment Length:</b> {top_hit.get('length', 0)} bp<br/>
        <b>E-value:</b> {top_hit.get('evalue', 0):.2e}<br/>
        <b>Genotype:</b> {genovar}
        """
        elements.append(Paragraph(blast_text, key_value_style))
    else:
        elements.append(Paragraph("No BLAST hits found", key_value_style))

    elements.append(Spacer(1, 0.1*inch))

    # MLST
    elements.append(Paragraph("MLST", subsection_style))
    mlst_text = f"""
    <b>Chlamydiales:</b> {data.get('mlst_chlamydiales', 'N/A')}<br/>
    <b>C. trachomatis:</b> {data.get('mlst_ct', 'N/A')}
    """
    elements.append(Paragraph(mlst_text, key_value_style))
    elements.append(Spacer(1, 0.1*inch))

    # Ct Strain Typing
    if ct:
        elements.append(Paragraph("Ct Strain Typing", subsection_style))
        ct_text = f"""
        <b>Assigned Strain:</b> {ct.get('type', 'N/A')}<br/>
        <b>Mash Best Match:</b> {ct.get('mash_best', 'N/A')}<br/>
        <b>FastANI Best Match:</b> {ct.get('fastani_best', 'N/A')}
        """
        elements.append(Paragraph(ct_text, key_value_style))

        # Interpretation
        elements.append(Spacer(1, 0.1*inch))
        ct_interp_style = ParagraphStyle(
            'CtInterpretation',
            parent=styles['Normal'],
            fontSize=9,
            textColor=colors.HexColor("#555555"),
            leftIndent=20
        )
        ct_interp_text = """
        <b>Method:</b> The sample genome was compared against 26 reference <i>C. trachomatis</i> genomes
        using Mash (k-mer distance) and FastANI (average nucleotide identity). The assigned strain
        represents the closest matching reference. FastANI values &gt;99.5% indicate very high similarity.<br/><br/>
        <b>Tip:</b> Users can add custom reference genomes to <font face="Courier">resources/references/ct/individual/</font>
        to improve strain typing for specific research needs.
        """
        elements.append(Paragraph(ct_interp_text, ct_interp_style))

    elements.append(Spacer(1, 0.2*inch))

    # =========================================================================
    # SECTION 3: PLASMID ANALYSIS
    # =========================================================================
    elements.append(Paragraph("3. Plasmid Analysis", section_style))

    elements.append(Paragraph("Dedicated Plasmid Assembly", subsection_style))

    # Status with color
    status = plasmid.get('status', 'Unknown')
    if status == "DETECTED":
        status_display = '<font color="#27AE60"><b>Ct PLASMID DETECTED</b></font>'
    elif status == "NOT_DETECTED":
        status_display = '<font color="#E74C3C"><b>NO PLASMID DETECTED</b></font>'
    else:
        status_display = f'<font color="#F39C12"><b>{status}</b></font>'

    plasmid_text = f"""
    <b>Detection Status:</b> {status_display}<br/>
    <b>Length:</b> {plasmid.get('length', 'N/A')} bp<br/>
    <b>Circular:</b> {plasmid.get('circular', 'Unknown')}<br/>
    <b>Size Deviation:</b> {plasmid.get('deviation', 'N/A')}
    """
    elements.append(Paragraph(plasmid_text, key_value_style))

    # Plasmid BLAST
    plasmid_blast = data.get("plasmid_blast", [])
    if plasmid_blast:
        elements.append(Spacer(1, 0.1*inch))
        elements.append(Paragraph("Plasmid BLAST Results", subsection_style))

        pblast_data = [["Subject", "Identity", "Length", "E-value"]]
        for hit in plasmid_blast[:5]:
            pblast_data.append([
                hit.get("subject", "")[:30],
                f"{hit.get('identity', 0):.1f}%",
                str(hit.get("length", 0)),
                f"{hit.get('evalue', 0):.1e}"
            ])

        pblast_table = Table(pblast_data, colWidths=[2.5*inch, 1*inch, 1*inch, 1.2*inch])
        pblast_table.setStyle(TableStyle([
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor("#27AE60")),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.white),
            ('FONTSIZE', (0, 0), (-1, -1), 8),
            ('BOX', (0, 0), (-1, -1), 1, colors.HexColor("#BDC3C7")),
            ('INNERGRID', (0, 0), (-1, -1), 0.5, colors.HexColor("#BDC3C7")),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 3),
            ('TOPPADDING', (0, 0), (-1, -1), 3),
        ]))
        elements.append(pblast_table)

    # Plasmid MLST
    plasmid_mlst = data.get("plasmid_mlst", "")
    if plasmid_mlst and plasmid_mlst != "Not available":
        elements.append(Spacer(1, 0.1*inch))
        elements.append(Paragraph(f"<b>Plasmid MLST:</b> {plasmid_mlst}", key_value_style))

    elements.append(Spacer(1, 0.2*inch))

    # =========================================================================
    # SECTION 4: GENOME ANNOTATION (if available)
    # =========================================================================
    annotation = data.get("annotation")
    if annotation:
        elements.append(Paragraph("4. Genome Annotation (Bakta)", section_style))

        # Calculate annotation efficiency metrics
        cds_count = int(annotation.get('cds', 0))
        hypothetical_count = int(annotation.get('hypotheticals', 0))
        annotated_count = cds_count - hypothetical_count if cds_count > hypothetical_count else cds_count
        annotation_rate = (annotated_count / cds_count * 100) if cds_count > 0 else 0

        annotation_data = [
            ["Feature", "Count"],
            ["Coding Sequences (CDS)", annotation.get('cds', 'N/A')],
            ["tRNAs", annotation.get('trnas', 'N/A')],
            ["rRNAs", annotation.get('rrnas', 'N/A')],
            ["ncRNAs", annotation.get('ncrnas', 'N/A')],
            ["Hypothetical Proteins", annotation.get('hypotheticals', 'N/A')],
            ["Pseudogenes", annotation.get('pseudogenes', 'N/A')],
        ]

        annotation_table = Table(annotation_data, colWidths=[3*inch, 2*inch])
        annotation_table.setStyle(TableStyle([
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor("#16A085")),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.white),
            ('FONTSIZE', (0, 0), (-1, -1), 9),
            ('BOX', (0, 0), (-1, -1), 1, colors.HexColor("#BDC3C7")),
            ('INNERGRID', (0, 0), (-1, -1), 0.5, colors.HexColor("#BDC3C7")),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 4),
            ('TOPPADDING', (0, 0), (-1, -1), 4),
            ('ALIGN', (1, 1), (1, -1), 'RIGHT'),
        ]))
        elements.append(annotation_table)

        # Annotation efficiency summary
        elements.append(Spacer(1, 0.1*inch))
        efficiency_text = f"""
        <b>Annotation Efficiency:</b> {annotation_rate:.1f}% of CDS have functional annotation<br/>
        <b>Functional CDS:</b> {annotated_count:,} / {cds_count:,}
        """
        elements.append(Paragraph(efficiency_text, key_value_style))

        elements.append(Spacer(1, 0.2*inch))

    # =========================================================================
    # SECTION 5: QUALITY CONTROL
    # =========================================================================
    section_num = "5" if annotation else "4"
    elements.append(Paragraph(f"{section_num}. Quality Control", section_style))

    # Read statistics
    fastp = data.get("fastp", {})
    if fastp:
        elements.append(Paragraph("Read Statistics", subsection_style))

        qc_data = [
            ["Metric", "Value"],
            ["Raw Reads", f"{fastp.get('raw_reads', 0):,}"],
            ["Filtered Reads", f"{fastp.get('filtered_reads', 0):,}"],
            ["Retention Rate", f"{100*fastp.get('filtered_reads', 0)/max(fastp.get('raw_reads', 1), 1):.1f}%"],
            ["Q30 Rate", f"{100*fastp.get('q30_rate', 0):.1f}%"],
            ["GC Content", f"{100*fastp.get('gc_content', 0):.1f}%"],
        ]

        qc_table = Table(qc_data, colWidths=[2.5*inch, 2.5*inch])
        qc_table.setStyle(TableStyle([
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor("#9B59B6")),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.white),
            ('FONTSIZE', (0, 0), (-1, -1), 9),
            ('BOX', (0, 0), (-1, -1), 1, colors.HexColor("#BDC3C7")),
            ('INNERGRID', (0, 0), (-1, -1), 0.5, colors.HexColor("#BDC3C7")),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 4),
            ('TOPPADDING', (0, 0), (-1, -1), 4),
            ('ALIGN', (1, 1), (1, -1), 'RIGHT'),
        ]))
        elements.append(qc_table)

    # Coverage
    coverage = data.get("coverage", {})
    if coverage or quast:
        elements.append(Spacer(1, 0.1*inch))
        elements.append(Paragraph("Coverage Statistics", subsection_style))

        # Get N50 from QUAST data (better metric for scaffolded single-chromosome assembly)
        n50_size = quast.get('n50', 0) if quast else 0

        # Calculate % genome coverage relative to standard Ct genome (1.043 Mb)
        CT_REFERENCE_SIZE = 1_043_000  # 1.043 Mb standard Ct genome size
        CT_MAX_EXPECTED = 1_060_000    # 1.06 Mb - flag assemblies larger than this
        genome_coverage_pct = (n50_size / CT_REFERENCE_SIZE * 100) if n50_size > 0 else 0

        cov_text = f"""
        <b>Mean Depth:</b> {coverage.get('mean_depth', 0):.1f}x<br/>
        <b>Contiguous Assembly (N50):</b> {n50_size:,} bp<br/>
        <b>Genome Coverage:</b> {genome_coverage_pct:.1f}%*
        """
        elements.append(Paragraph(cov_text, key_value_style))

        # QC flag for assemblies larger than expected (>1.06 Mb)
        if n50_size > CT_MAX_EXPECTED:
            warning_style = ParagraphStyle(
                'Warning',
                parent=normal_style,
                fontSize=9,
                textColor=colors.HexColor("#E74C3C"),
                leftIndent=20,
                spaceBefore=4
            )
            warning_text = f"""
            <b>⚠ QC FLAG:</b> Assembly N50 ({n50_size:,} bp) exceeds expected Ct genome size (>1.06 Mb).<br/>
            Review for: contamination, uncollapsed duplications, or scaffolding artifacts.
            """
            elements.append(Paragraph(warning_text, warning_style))

        # Footnote for genome coverage assumption
        footnote_style = ParagraphStyle(
            'Footnote',
            parent=normal_style,
            fontSize=8,
            textColor=colors.HexColor("#7F8C8D"),
            leftIndent=20
        )
        elements.append(Paragraph("*N50 relative to standard C. trachomatis genome size (1.043 Mb)", footnote_style))

    # Kraken
    kraken = data.get("kraken", [])
    if kraken:
        elements.append(Spacer(1, 0.1*inch))
        elements.append(Paragraph("Kraken2 Classification (Top Hits)", subsection_style))

        kraken_data = [["Taxon", "Percentage"]]
        for hit in kraken[:5]:
            kraken_data.append([hit.get("name", "")[:40], f"{hit.get('percent', 0):.1f}%"])

        kraken_table = Table(kraken_data, colWidths=[4*inch, 1.5*inch])
        kraken_table.setStyle(TableStyle([
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor("#E67E22")),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.white),
            ('FONTSIZE', (0, 0), (-1, -1), 8),
            ('BOX', (0, 0), (-1, -1), 1, colors.HexColor("#BDC3C7")),
            ('INNERGRID', (0, 0), (-1, -1), 0.5, colors.HexColor("#BDC3C7")),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 3),
            ('TOPPADDING', (0, 0), (-1, -1), 3),
            ('ALIGN', (1, 1), (1, -1), 'RIGHT'),
        ]))
        elements.append(kraken_table)

    # =========================================================================
    # SECTION: MIXED STRAIN ANALYSIS
    # =========================================================================
    mixed = data.get("mixed_detection")
    if mixed:
        # Determine section number
        mixed_section_num = "6" if annotation else "5"
        elements.append(Spacer(1, 0.2*inch))
        elements.append(Paragraph(f"{mixed_section_num}. Mixed Strain Analysis", section_style))

        het_sites = mixed.get('het_sites', 0)
        mean_af = mixed.get('mean_af', 'N/A')
        status = mixed.get('status', 'Unknown')
        main_contig_name = mixed.get('main_contig_name', 'N/A')
        main_contig_size = mixed.get('main_contig_size', 0)
        position_range = mixed.get('position_range', 0)
        distribution = mixed.get('distribution', 'N/A')

        # Format sizes
        if position_range > 1000:
            range_str = f"{position_range/1000:.1f} kb"
        else:
            range_str = f"{position_range} bp"

        if main_contig_size > 1000:
            contig_size_str = f"{main_contig_size/1000:.1f} kb"
        else:
            contig_size_str = f"{main_contig_size} bp"

        # Summary table with key metrics
        mixed_data = [
            ["Metric", "Value", "Interpretation"],
            ["Heterozygous Sites", str(het_sites), "Sites passing all quality filters"],
            ["Mean Allele Frequency", str(mean_af), "Expected ~0.5 for equal strain mixture"],
            ["Main Chromosome", contig_size_str, "Largest contig analyzed (excludes plasmid/fragments)"],
            ["Position Spread", range_str, "Genomic range of variants on main chromosome"],
        ]

        mixed_table = Table(mixed_data, colWidths=[1.8*inch, 1.2*inch, 3*inch])
        mixed_table.setStyle(TableStyle([
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor("#8E44AD")),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.white),
            ('FONTSIZE', (0, 0), (-1, -1), 9),
            ('BOX', (0, 0), (-1, -1), 1, colors.HexColor("#BDC3C7")),
            ('INNERGRID', (0, 0), (-1, -1), 0.5, colors.HexColor("#BDC3C7")),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 4),
            ('TOPPADDING', (0, 0), (-1, -1), 4),
            ('ALIGN', (1, 1), (1, -1), 'CENTER'),
        ]))
        elements.append(mixed_table)
        elements.append(Spacer(1, 0.1*inch))

        # Status interpretation
        if status == "MIXED_WARNING":
            status_color = "#E74C3C"
            dist_text = ""
            if distribution == "GENOME-WIDE":
                dist_text = "<br/><br/><b>Distribution:</b> <font color='#27AE60'>GENOME-WIDE</font> - Variants are spread across the genome, consistent with true mixed infection."
            elif distribution == "LOCALIZED":
                dist_text = "<br/><br/><b>Distribution:</b> <font color='#F39C12'>LOCALIZED</font> - Variants are clustered in one region, could indicate recombination or localized artifact."

            status_text = f"""
            <b>Status:</b> <font color="{status_color}">POSSIBLE MIXED INFECTION</font><br/><br/>
            <b>Interpretation:</b> This sample shows {het_sites} genomic positions where sequencing reads
            support two different alleles at intermediate frequencies. Since <i>C. trachomatis</i> is haploid,
            this pattern suggests the presence of multiple distinct strains in the sample.{dist_text}<br/><br/>
            <b>Recommendations:</b><br/>
            • Interpret typing results (genotype, MLST, Ct strain) with caution - they may represent a consensus<br/>
            • Consider that plasmid detection may be affected by strain mixture<br/>
            • Review clinical context for possible co-infection or reinfection
            """
        elif status == "SINGLE" or status == "LOW_HET":
            status_color = "#27AE60"
            status_text = f"""
            <b>Status:</b> <font color="{status_color}">SINGLE STRAIN (LOW HETEROZYGOSITY)</font><br/><br/>
            <b>Interpretation:</b> This sample shows minimal heterozygosity ({het_sites} sites),
            consistent with a single <i>C. trachomatis</i> strain. Typing results can be interpreted
            with normal confidence.
            """
        else:
            status_color = "#F39C12"
            status_text = f"""
            <b>Status:</b> <font color="{status_color}">INCONCLUSIVE</font><br/><br/>
            <b>Interpretation:</b> Mixed strain analysis could not determine a clear result.
            """

        elements.append(Paragraph(status_text, key_value_style))

        # Methodology note
        elements.append(Spacer(1, 0.1*inch))
        method_note_style = ParagraphStyle(
            'MethodNote',
            parent=normal_style,
            fontSize=8,
            textColor=colors.HexColor("#7F8C8D"),
            leftIndent=20
        )
        method_text = """
        <b>Method:</b> Reads were mapped to the sample's main chromosome (largest contig) and variants called.
        Heterozygous-like sites are positions where 15-85% of reads support an alternate allele.
        <b>Quality filters:</b> Base quality ≥20, minimum 2 reads on each strand (reduces sequencing errors
        and strand bias artifacts). Analysis focuses on main chromosome only (excludes plasmid and small fragments).
        In a haploid organism, sites passing these filters indicate multiple genomes.
        """
        elements.append(Paragraph(method_text, method_note_style))

    # =========================================================================
    # FOOTER
    # =========================================================================
    elements.append(Spacer(1, 0.3*inch))
    elements.append(HRFlowable(width="100%", thickness=1, color=colors.HexColor("#BDC3C7")))
    elements.append(Spacer(1, 0.1*inch))

    footer_text = f"""
    <font size=8 color="#7F8C8D">
    Generated by CtGAP (Chlamydia trachomatis Genome Analysis Pipeline)<br/>
    Report created: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}
    </font>
    """
    elements.append(Paragraph(footer_text, ParagraphStyle('Footer', alignment=TA_CENTER)))

    # Build PDF
    doc.build(elements)


def main():
    # Get inputs from snakemake
    sample = snakemake.params.sample
    ct_typing_enabled = snakemake.params.ct_typing_enabled

    # Collect all data
    data = {}

    # Parse each input
    # First get assembly method to determine correct QUAST source
    data["assembly_selection"] = parse_assembly_selection(
        snakemake.input.assembly_report,
        snakemake.input.assembly_method
    )

    # Use QUAST from the winning assembly method (not always denovo)
    assembly_method = data["assembly_selection"].get("method", "denovo")
    outdir = Path(snakemake.params.outdir)
    if assembly_method == "reference-denovo":
        quast_path = outdir / sample / "ref-denovo" / "quast" / "report.tsv"
    else:
        quast_path = outdir / sample / "denovo" / "assembly_statistics" / "report.tsv"
    data["quast"] = parse_quast_report(str(quast_path))

    # Calculate gaps (N's) in the main contig from the best assembly
    assembly_file = str(snakemake.input.best_assembly)
    data["quast"]["main_contig_gaps"] = count_gaps_in_main_contig(assembly_file)

    data["blast_ompa"] = parse_blast_results(snakemake.input.blast_ompa)
    data["blast_secondary"] = parse_blast_results(snakemake.input.blast_secondary)
    data["mlst_chlamydiales"] = parse_mlst_results(snakemake.input.mlst_chlamydiales)
    data["mlst_ct"] = parse_mlst_results(snakemake.input.mlst_ct)
    # assembly_selection already parsed above for QUAST path determination
    data["plasmid_stats"] = parse_plasmid_analysis(snakemake.input.plasmid_analysis)
    data["plasmid_blast"] = parse_blast_results(snakemake.input.plasmid_blast)
    data["plasmid_mlst"] = parse_mlst_results(snakemake.input.plasmid_mlst)
    data["fastp"] = parse_fastp_json(snakemake.input.fastp_json)
    data["kraken"] = parse_kraken_report(snakemake.input.kraken_report)
    data["coverage"] = parse_coverage(snakemake.input.coverage)

    # CT typing (if enabled)
    if ct_typing_enabled and snakemake.input.ct_typing:
        data["ct_typing"] = parse_ct_typing(snakemake.input.ct_typing)
    else:
        data["ct_typing"] = None

    # Annotation (if available - check for file existence)
    annotation_enabled = snakemake.params.annotation_enabled
    # outdir already defined above for QUAST path
    if annotation_enabled:
        annotation_txt = outdir / sample / "annotation" / f"{sample}.txt"
        data["annotation"] = parse_bakta_annotation(str(annotation_txt))
    else:
        data["annotation"] = None

    # Mixed strain detection (if available)
    mixed_detection_enabled = snakemake.params.mixed_detection_enabled
    if mixed_detection_enabled:
        mixed_summary = outdir / sample / "mixed_detection" / f"{sample}.mixed_summary.txt"
        data["mixed_detection"] = parse_mixed_detection(str(mixed_summary))
    else:
        data["mixed_detection"] = None

    # Create PDF
    create_pdf_report(sample, snakemake.output.pdf, data)

    # Touch status file
    Path(snakemake.output.status).touch()


if __name__ == "__main__":
    main()
