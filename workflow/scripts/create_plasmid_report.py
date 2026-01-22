#!/usr/bin/env python3
"""
Create plasmid analysis PDF report for CtGAP.
Uses reference-guided read extraction + Unicycler assembly approach.
"""

import os
import sys
from pathlib import Path
from datetime import datetime

from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, HRFlowable
)


def parse_analysis_file(filepath):
    """Parse the plasmid analysis text file."""
    result = {
        "status": "UNKNOWN",
        "length": "",
        "circular": "",
        "deviation": ""
    }

    if not os.path.exists(filepath):
        return result

    with open(filepath) as f:
        content = f.read()

    if "CT PLASMID DETECTED" in content:
        result["status"] = "DETECTED"
    elif "NO PLASMID DETECTED" in content:
        result["status"] = "NOT_DETECTED"
    elif "NO VALID CT PLASMID" in content:
        result["status"] = "NO_VALID_PLASMID"

    for line in content.split("\n"):
        if "Length:" in line and "bp" in line:
            result["length"] = line.split(":")[1].strip()
        if "Circular:" in line:
            result["circular"] = line.split(":")[1].strip()
        if "Size deviation:" in line:
            result["deviation"] = line.split(":")[1].strip()

    return result


def parse_blast_results(filepath):
    """Parse BLAST tabular output."""
    results = []

    if not os.path.exists(filepath) or os.path.getsize(filepath) == 0:
        return results

    with open(filepath) as f:
        for line in f:
            if line.startswith("query"):
                continue
            fields = line.strip().split("\t")
            if len(fields) >= 12:
                results.append({
                    "subject": fields[1],
                    "identity": float(fields[2]),
                    "length": int(fields[3]),
                    "evalue": float(fields[10]),
                    "bitscore": float(fields[11])
                })

    return sorted(results, key=lambda x: x["bitscore"], reverse=True)


def parse_mapping_stats(filepath):
    """Parse samtools flagstat output."""
    result = {"mapped_reads": 0, "mapped_percent": "0%"}

    if not os.path.exists(filepath):
        return result

    with open(filepath) as f:
        for line in f:
            if "mapped (" in line and "primary" not in line:
                parts = line.split()
                if parts:
                    result["mapped_reads"] = int(parts[0])
                    # Extract percentage
                    if "(" in line and "%" in line:
                        pct = line.split("(")[1].split("%")[0]
                        result["mapped_percent"] = f"{pct}%"

    return result


def parse_mlst_results(filepath):
    """Parse MLST output."""
    if not os.path.exists(filepath) or os.path.getsize(filepath) == 0:
        return "Not available"

    with open(filepath) as f:
        lines = [l.strip() for l in f.readlines() if l.strip()]

    if len(lines) > 1:
        # Skip header, get data line
        fields = lines[1].split("\t")
        if len(fields) >= 2:
            return f"ST{fields[1]}" if fields[1] != "NA" else "No match"

    return "No match"


def create_report(sample, output_pdf, analysis, blast, mlst, mapping):
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

    elements = []

    # Title
    elements.append(Paragraph("CT Plasmid Analysis Report", title_style))
    elements.append(Spacer(1, 0.1*inch))

    # Header info
    header_data = [
        ["Sample ID:", sample],
        ["Report Date:", datetime.now().strftime("%Y-%m-%d %H:%M")],
        ["Method:", "Reference-guided read extraction + Unicycler"],
    ]
    header_table = Table(header_data, colWidths=[1.5*inch, 4*inch])
    header_table.setStyle(TableStyle([
        ('FONTNAME', (0, 0), (0, -1), 'Helvetica-Bold'),
        ('FONTSIZE', (0, 0), (-1, -1), 10),
        ('BOTTOMPADDING', (0, 0), (-1, -1), 4),
    ]))
    elements.append(header_table)
    elements.append(Spacer(1, 0.2*inch))

    elements.append(HRFlowable(width="100%", thickness=2, color=colors.HexColor("#3498DB")))
    elements.append(Spacer(1, 0.2*inch))

    # Main Result Box
    elements.append(Paragraph("Plasmid Detection Result", section_style))

    # Status color
    if analysis["status"] == "DETECTED":
        status_color = colors.HexColor("#27AE60")  # Green
        status_text = "PLASMID DETECTED"
    elif analysis["status"] == "NOT_DETECTED":
        status_color = colors.HexColor("#E74C3C")  # Red
        status_text = "NO PLASMID DETECTED"
    else:
        status_color = colors.HexColor("#F39C12")  # Orange
        status_text = "UNCERTAIN"

    result_data = [
        ["Status:", status_text],
        ["Length:", analysis.get("length", "N/A")],
        ["Circular:", analysis.get("circular", "Unknown")],
        ["Size Deviation:", analysis.get("deviation", "N/A")],
    ]

    result_table = Table(result_data, colWidths=[1.5*inch, 4*inch])
    result_table.setStyle(TableStyle([
        ('FONTNAME', (0, 0), (0, -1), 'Helvetica-Bold'),
        ('FONTSIZE', (0, 0), (-1, -1), 11),
        ('BACKGROUND', (1, 0), (1, 0), status_color),
        ('TEXTCOLOR', (1, 0), (1, 0), colors.white),
        ('BOTTOMPADDING', (0, 0), (-1, -1), 6),
        ('TOPPADDING', (0, 0), (-1, -1), 6),
        ('BOX', (0, 0), (-1, -1), 1, colors.HexColor("#BDC3C7")),
    ]))
    elements.append(result_table)
    elements.append(Spacer(1, 0.2*inch))

    # Mapping Statistics
    elements.append(Paragraph("Read Mapping to Reference Plasmid", section_style))

    mapping_data = [
        ["Metric", "Value"],
        ["Reads mapped to pCT reference", f"{mapping.get('mapped_reads', 0):,}"],
        ["Mapping rate", mapping.get("mapped_percent", "0%")],
    ]

    mapping_table = Table(mapping_data, colWidths=[3*inch, 2*inch])
    mapping_table.setStyle(TableStyle([
        ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
        ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor("#9B59B6")),
        ('TEXTCOLOR', (0, 0), (-1, 0), colors.white),
        ('FONTSIZE', (0, 0), (-1, -1), 10),
        ('BOX', (0, 0), (-1, -1), 1, colors.HexColor("#BDC3C7")),
        ('INNERGRID', (0, 0), (-1, -1), 0.5, colors.HexColor("#BDC3C7")),
        ('BOTTOMPADDING', (0, 0), (-1, -1), 4),
        ('TOPPADDING', (0, 0), (-1, -1), 4),
    ]))
    elements.append(mapping_table)
    elements.append(Spacer(1, 0.2*inch))

    # BLAST Results
    elements.append(Paragraph("BLAST Results (vs 20 CT Reference Plasmids)", section_style))

    if blast:
        blast_data = [["Reference Plasmid", "Identity (%)", "Alignment", "E-value"]]
        for hit in blast[:5]:
            blast_data.append([
                hit["subject"][:30],
                f"{hit['identity']:.1f}",
                f"{hit['length']} bp",
                f"{hit['evalue']:.1e}"
            ])

        blast_table = Table(blast_data, colWidths=[2.5*inch, 1*inch, 1*inch, 1*inch])
        blast_table.setStyle(TableStyle([
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor("#27AE60")),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.white),
            ('FONTSIZE', (0, 0), (-1, -1), 9),
            ('BOX', (0, 0), (-1, -1), 1, colors.HexColor("#BDC3C7")),
            ('INNERGRID', (0, 0), (-1, -1), 0.5, colors.HexColor("#BDC3C7")),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 3),
            ('TOPPADDING', (0, 0), (-1, -1), 3),
        ]))
        elements.append(blast_table)
    else:
        elements.append(Paragraph("No BLAST hits found.", styles['Normal']))

    elements.append(Spacer(1, 0.2*inch))

    # MLST
    elements.append(Paragraph("Plasmid MLST", section_style))
    elements.append(Paragraph(f"<b>MLST Result:</b> {mlst}", styles['Normal']))

    elements.append(Spacer(1, 0.3*inch))

    # Methodology note
    elements.append(HRFlowable(width="100%", thickness=1, color=colors.HexColor("#BDC3C7")))
    elements.append(Spacer(1, 0.1*inch))

    method_text = """
    <font size=8 color="#7F8C8D">
    <b>Methodology:</b> Reads were mapped to the reference CT plasmid (pCT, D strain, ~7,500 bp).
    Mapped read pairs were extracted and assembled using Unicycler. The assembly was verified
    by size (expected 6,000-10,000 bp) and BLAST against 20 CT reference plasmids.<br/><br/>
    Generated by CtGAP (Chlamydia trachomatis Genome Analysis Pipeline)
    </font>
    """
    elements.append(Paragraph(method_text, styles['Normal']))

    # Build PDF
    doc.build(elements)


def main():
    sample = snakemake.wildcards.sample

    analysis = parse_analysis_file(snakemake.input.analysis)
    blast = parse_blast_results(snakemake.input.blast)
    mlst = parse_mlst_results(snakemake.input.mlst)
    mapping = parse_mapping_stats(snakemake.input.mapping_stats)

    create_report(
        sample,
        snakemake.output.pdf,
        analysis,
        blast,
        mlst,
        mapping
    )


if __name__ == "__main__":
    main()
