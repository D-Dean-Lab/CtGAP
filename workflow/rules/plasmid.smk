# =============================================================================
# PLASMID ANALYSIS MODULE 
# =============================================================================
# Uses Unicycler for bacterial genome + plasmid assembly and separation
# =============================================================================

# -----------------------------------------------------------------------------
# 1. REFERENCE PLASMID EXTRACTION
# -----------------------------------------------------------------------------
rule extract_reference_plasmid:
    input:
        plasmids = "resources/references/ct/20_Ct_plasmids.fasta"
    output:
        reference = "resources/references/ct/reference_plasmid.fasta"
    log: "logs/extract_reference_plasmid.log"
    shell: r"""
        echo "Extracting plasmid D from $(basename {input.plasmids})" > {log}

        # Extract the D plasmid sequence
        awk '
        /^>.*D/ {{
            found = 1
            print $0
            next
        }}
        found && /^>/ {{
            found = 0
            exit
        }}
        found {{
            print $0
        }}
        ' {input.plasmids} > {output.reference}

        # Check if D plasmid was found
        if [ ! -s {output.reference} ]; then
            echo "ERROR: Plasmid D not found in {input.plasmids}" >> {log}
            echo "Available sequences:" >> {log}
            grep '^>' {input.plasmids} >> {log}
            exit 1
        fi

        echo "Plasmid D extracted successfully" >> {log}
        echo "Reference length: $(grep -v '^>' {output.reference} | tr -d '\n' | wc -c) bp" >> {log}
        echo "Reference name: $(grep '^>' {output.reference})" >> {log}
    """

# -----------------------------------------------------------------------------
# 2. UNICYCLER ASSEMBLY 
# -----------------------------------------------------------------------------
rule unicycler_assembly:
    input:
        r1 = rules.extract_chlamydiales_reads.output.r1,
        r2 = rules.extract_chlamydiales_reads.output.r2,
    output:
        assembly = OUTDIR / "{sample}" / "denovo" / "unicycler" / "assembly.fasta",
        graph = OUTDIR / "{sample}" / "denovo" / "unicycler" / "assembly.gfa", 
        log_file = OUTDIR / "{sample}" / "denovo" / "unicycler" / "unicycler.log",
        status = OUTDIR / "status" / "denovo.unicycler.{sample}.txt"
    params:
        outdir = lambda wildcards: str(OUTDIR / wildcards.sample / "denovo" / "unicycler")
    threads: config['threads']['spades']  # Use existing spades thread config
    conda: "../envs/unicycler.yaml"
    log: OUTDIR / "{sample}" / "log" / "denovo.unicycler.{sample}.log"
    benchmark: OUTDIR / "{sample}" / "benchmark" / "unicycler.{sample}.txt"
    shell: r"""
        # Remove existing output directory
        rm -rf {params.outdir}
        mkdir -p {params.outdir}
        
        # Run Unicycler
        unicycler \
            -1 {input.r1} \
            -2 {input.r2} \
            -o {params.outdir} \
            -t {threads} \
            --verbosity 2 \
            > {log} 2>&1
        
        # Check if assembly was successful
        if [ -f "{output.assembly}" ] && [ -s "{output.assembly}" ]; then
            echo "UNICYCLER_SUCCESS" > {params.outdir}/UNICYCLER_SUCCESS.txt
            echo "Unicycler assembly completed successfully" >> {log}
            
            # Log assembly statistics
            echo "Assembly statistics:" >> {log}
            echo "Total contigs: $(grep -c '^>' {output.assembly})" >> {log}
            echo "Contig sizes:" >> {log}
            awk '/^>/ {{if(seq) print length(seq) " bp"; seq=""; next}} {{seq=seq$0}} END {{if(seq) print length(seq) " bp"}}' {output.assembly} >> {log}
        else
            echo "ERROR: Unicycler assembly failed" >> {log}
            echo "UNICYCLER_FAILED" > {params.outdir}/UNICYCLER_FAILED.txt
            touch {output.assembly} {output.graph} {output.log_file}
        fi
        
        touch {output.status}
    """

# -----------------------------------------------------------------------------
# 3. EXTRACT PLASMIDS FROM UNICYCLER OUTPUT (with copy number calculation)
# -----------------------------------------------------------------------------
rule extract_unicycler_plasmids:
    input:
        assembly = rules.unicycler_assembly.output.assembly,
        graph = rules.unicycler_assembly.output.graph
    output:
        plasmids = OUTDIR / "{sample}" / "denovo" / "unicycler" / "plasmids_only.fasta",
        chromosome = OUTDIR / "{sample}" / "denovo" / "unicycler" / "chromosome.fasta",
        stats = OUTDIR / "{sample}" / "denovo" / "unicycler" / "assembly_stats.txt",
        status = OUTDIR / "status" / "denovo.unicycler_plasmids.{sample}.txt"
    shell: r"""
        # Create stats file
        echo "Unicycler Assembly Analysis for {wildcards.sample}" > {output.stats}
        echo "Generated: $(date)" >> {output.stats}
        echo "==========================================" >> {output.stats}
        echo "" >> {output.stats}
        
        if [ -s {input.assembly} ]; then
            # Initialize temp files
            temp_dir=$(mktemp -d)
            plasmid_seqs="$temp_dir/plasmids.fa"
            chr_seqs="$temp_dir/chromosome.fa"
            
            echo "Contig analysis:" >> {output.stats}
            echo "---------------" >> {output.stats}
            
            # Use Python for more reliable parsing (simpler than complex AWK)
            python3 -c "
import sys
import re

plasmid_count = 0
chr_count = 0
chr_depth = 1.0

# Read the assembly and process each contig
with open('{input.assembly}', 'r') as f:
    header = ''
    seq = ''
    
    for line in f:
        line = line.strip()
        if line.startswith('>'):
            # Process previous contig if it exists
            if header and seq:
                length = len(seq)
                contig_name = header[1:].split()[0]  # Remove > and get first part
                
                # Extract depth information using regex
                depth_match = re.search(r'depth=([0-9]+\.?[0-9]*)x', header)
                current_depth = float(depth_match.group(1)) if depth_match else 1.0
                
                # Classification logic
                if length >= 1000 and length <= 500000:
                    copy_number = current_depth / chr_depth
                    if copy_number >= 50:
                        copy_desc = f' (~{{int(copy_number + 0.5)}} copies - HIGH COPY)'
                    elif copy_number >= 10:
                        copy_desc = f' (~{{int(copy_number + 0.5)}} copies - MEDIUM COPY)'
                    else:
                        copy_desc = f' (~{{int(copy_number + 0.5)}} copies - LOW COPY)'
                    
                    if length <= 50000:
                        classification = '[PLASMID - high confidence]' + copy_desc
                        with open('$plasmid_seqs', 'a') as pf:
                            pf.write(header + '\\n' + seq + '\\n')
                        plasmid_count += 1
                    else:
                        classification = '[POSSIBLE_PLASMID - check manually]' + copy_desc
                        with open('$plasmid_seqs', 'a') as pf:
                            pf.write(header + '\\n' + seq + '\\n')
                        plasmid_count += 1
                        
                elif length >= 500000:
                    chr_depth = current_depth  # Update chromosome depth
                    classification = '[CHROMOSOME]'
                    with open('$chr_seqs', 'a') as cf:
                        cf.write(header + '\\n' + seq + '\\n')
                    chr_count += 1
                else:
                    classification = '[TOO_SMALL - likely artifact]'
                
                # Write to stats file
                with open('{output.stats}', 'a') as sf:
                    sf.write(f'  {{contig_name}}: {{length}} bp {{classification}}\\n')
            
            # Start new contig
            header = line
            seq = ''
        else:
            seq += line
    
    # Process last contig
    if header and seq:
        length = len(seq)
        contig_name = header[1:].split()[0]
        
        depth_match = re.search(r'depth=([0-9]+\.?[0-9]*)x', header)
        current_depth = float(depth_match.group(1)) if depth_match else 1.0
        
        if length >= 1000 and length <= 500000:
            copy_number = current_depth / chr_depth
            if copy_number >= 50:
                copy_desc = f' (~{{int(copy_number + 0.5)}} copies - HIGH COPY)'
            elif copy_number >= 10:
                copy_desc = f' (~{{int(copy_number + 0.5)}} copies - MEDIUM COPY)'
            else:
                copy_desc = f' (~{{int(copy_number + 0.5)}} copies - LOW COPY)'
            
            if length <= 50000:
                classification = '[PLASMID - high confidence]' + copy_desc
                with open('$plasmid_seqs', 'a') as pf:
                    pf.write(header + '\\n' + seq + '\\n')
                plasmid_count += 1
            else:
                classification = '[POSSIBLE_PLASMID - check manually]' + copy_desc
                with open('$plasmid_seqs', 'a') as pf:
                    pf.write(header + '\\n' + seq + '\\n')
                plasmid_count += 1
                
        elif length >= 500000:
            chr_depth = current_depth
            classification = '[CHROMOSOME]'
            with open('$chr_seqs', 'a') as cf:
                cf.write(header + '\\n' + seq + '\\n')
            chr_count += 1
        else:
            classification = '[TOO_SMALL - likely artifact]'
        
        with open('{output.stats}', 'a') as sf:
            sf.write(f'  {{contig_name}}: {{length}} bp {{classification}}\\n')

# Write summary
with open('{output.stats}', 'a') as sf:
    sf.write('\\n')
    sf.write('Summary:\\n')
    sf.write(f'  Plasmids found: {{plasmid_count}}\\n')
    sf.write(f'  Chromosomes found: {{chr_count}}\\n')
    
    if plasmid_count > 0:
        sf.write('\\nSUCCESS: Plasmids found!\\n')
    else:
        sf.write('\\nWARNING: No plasmids detected in this sample\\n')
"
            
            # Move results to final outputs
            if [ -s "$plasmid_seqs" ]; then
                cp "$plasmid_seqs" {output.plasmids}
                echo "Plasmid extraction completed!"
            else
                echo "No plasmids found - creating empty file" >> {output.stats}
                touch {output.plasmids}
            fi
            
            if [ -s "$chr_seqs" ]; then
                cp "$chr_seqs" {output.chromosome}
            else
                touch {output.chromosome}
            fi
            
            # Clean up temp files
            rm -rf "$temp_dir"
        else
            echo "ERROR: Input assembly file is empty" >> {output.stats}
            touch {output.plasmids}
            touch {output.chromosome}
        fi

        touch {output.status}
    """

# -----------------------------------------------------------------------------
# 4. EXTRACT MAIN PLASMID (largest plasmid for downstream analysis) 
# -----------------------------------------------------------------------------
rule extract_main_plasmid:
    input:
        plasmids = rules.extract_unicycler_plasmids.output.plasmids
    output:
        main_plasmid = OUTDIR / "{sample}" / "denovo" / "unicycler" / "{sample}.main_plasmid.fasta",
        status = OUTDIR / "status" / "denovo.extract_main_plasmid.{sample}.txt"
    shell: r"""
        if [ -s {input.plasmids} ]; then
            # Extract the largest plasmid (most likely to be the main one)
            awk '
            /^>/ {{
                if(seq && header && length(seq) > max_len) {{
                    max_len = length(seq)
                    max_header = header
                    max_seq = seq
                }}
                header = $0
                seq = ""
                next
            }}
            {{seq = seq $0}}
            END {{
                if(seq && header && length(seq) > max_len) {{
                    max_len = length(seq)
                    max_header = header  
                    max_seq = seq
                }}
                if(max_header) {{
                    print max_header
                    print max_seq
                }}
            }}
            ' {input.plasmids} > {output.main_plasmid}
        else
            echo "No plasmids available for extraction" 
            touch {output.main_plasmid}
        fi

        touch {output.status}
    """

# -----------------------------------------------------------------------------
# 5. CREATE SAMPLE PDF REPORT
# -----------------------------------------------------------------------------
rule create_plasmid_pdf_report:
    input:
        assembly_stats = rules.extract_unicycler_plasmids.output.stats,
        plasmids_fasta = rules.extract_unicycler_plasmids.output.plasmids,
        main_plasmid = rules.extract_main_plasmid.output.main_plasmid,
        blast_results = OUTDIR / "{sample}" / "denovo" / "blast" / "blast.plasmid.tab",
        mlst_results = OUTDIR / "{sample}" / "denovo" / "mlst" / "{sample}.genome.plasmid.mlst.txt"
    output:
        pdf_report = OUTDIR / "{sample}" / "reports" / "{sample}_plasmid_report.pdf"
    conda: "../envs/pdf.yaml"
    log: OUTDIR / "{sample}" / "log" / "plasmid_pdf_report.{sample}.log"
    shell: r"""
        python3 -c "
import os
import sys
from pathlib import Path
from datetime import datetime
from reportlab.lib.pagesizes import letter, A4
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, PageBreak, Table, TableStyle
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch
from reportlab.lib import colors
from reportlab.pdfgen import canvas

def parse_assembly_stats(stats_file):
    '''Parse the assembly stats file'''
    data = {{}}
    plasmids = []
    chromosomes = []
    
    try:
        with open(stats_file, 'r') as f:
            content = f.read()
            
        # Extract sample name and date
        lines = content.split('\n')
        for line in lines:
            if line.startswith('Unicycler Assembly Analysis for'):
                data['sample'] = line.split('for ')[-1].strip()
            elif line.startswith('Generated:'):
                data['date'] = line.split('Generated: ')[-1].strip()
            elif 'bp [PLASMID' in line:
                plasmids.append(line.strip())
            elif 'bp [CHROMOSOME]' in line:
                chromosomes.append(line.strip())
            elif 'Plasmids found:' in line:
                data['plasmid_count'] = line.split(': ')[-1].strip()
            elif 'Chromosomes found:' in line:
                data['chromosome_count'] = line.split(': ')[-1].strip()
                
        data['plasmids'] = plasmids
        data['chromosomes'] = chromosomes
        
    except Exception as e:
        print(f'Error parsing assembly stats: {{e}}', file=sys.stderr)
        data = {{'sample': '{wildcards.sample}', 'date': 'Unknown', 'plasmid_count': '0', 'chromosome_count': '0', 'plasmids': [], 'chromosomes': []}}
    
    return data

def parse_blast_results(blast_file):
    '''Parse BLAST results'''
    blast_data = []
    try:
        if os.path.exists(blast_file) and os.path.getsize(blast_file) > 0:
            # BLAST format: query subject pident length mismatch gapopen qstart qend sstart send evalue bitscore
            with open(blast_file, 'r') as f:
                for line in f:
                    if line.strip():
                        fields = line.strip().split('\t')
                        if len(fields) >= 12:
                            blast_data.append({{
                                'subject': fields[1],
                                'identity': float(fields[2]),
                                'length': int(fields[3]),
                                'evalue': float(fields[10]),
                                'bitscore': float(fields[11])
                            }})
    except Exception as e:
        print(f'Error parsing BLAST results: {{e}}', file=sys.stderr)
    
    return blast_data

def parse_mlst_results(mlst_file):
    '''Parse MLST results'''
    mlst_data = {{}}
    try:
        if os.path.exists(mlst_file) and os.path.getsize(mlst_file) > 0:
            with open(mlst_file, 'r') as f:
                content = f.read().strip()
                if content:
                    mlst_data['result'] = content
                else:
                    mlst_data['result'] = 'No MLST match found'
        else:
            mlst_data['result'] = 'No MLST data available'
    except Exception as e:
        print(f'Error parsing MLST results: {{e}}', file=sys.stderr)
        mlst_data = {{'result': 'Error reading MLST data'}}
    
    return mlst_data

def count_plasmid_sequences(fasta_file):
    '''Count sequences in plasmid FASTA'''
    count = 0
    total_length = 0
    try:
        if os.path.exists(fasta_file) and os.path.getsize(fasta_file) > 0:
            with open(fasta_file, 'r') as f:
                seq_len = 0
                for line in f:
                    if line.startswith('>'):
                        if seq_len > 0:
                            total_length += seq_len
                            seq_len = 0
                        count += 1
                    else:
                        seq_len += len(line.strip())
                # Add the last sequence
                if seq_len > 0:
                    total_length += seq_len
    except Exception as e:
        print(f'Error counting sequences: {{e}}', file=sys.stderr)
        
    return count, total_length

def create_pdf_report(output_path, assembly_data, blast_data, mlst_data, plasmid_count, plasmid_length):
    '''Create the main PDF report'''
    
    doc = SimpleDocTemplate(output_path, pagesize=A4,
                           rightMargin=72, leftMargin=72,
                           topMargin=72, bottomMargin=18)
    
    # Get styles
    styles = getSampleStyleSheet()
    title_style = ParagraphStyle(
        'CustomTitle',
        parent=styles['Heading1'],
        fontSize=20,
        spaceAfter=30,
        textColor=colors.darkblue
    )
    
    heading_style = ParagraphStyle(
        'CustomHeading',
        parent=styles['Heading2'], 
        fontSize=14,
        spaceAfter=12,
        spaceBefore=20,
        textColor=colors.darkred
    )
    
    # Build the story
    story = []
    
    # Title page
    title = Paragraph(f'Plasmid Analysis Report', title_style)
    story.append(title)
    
    subtitle = Paragraph(f'Sample: {{assembly_data.get(\"sample\", \"Unknown\")}}', styles['Heading2'])
    story.append(subtitle)
    story.append(Spacer(1, 20))
    
    # Report info
    info_data = [
        ['Generated:', assembly_data.get('date', 'Unknown')],
        ['Pipeline:', 'CtGAP 2.0'],
        ['Analysis Type:', 'Plasmid Assembly & Characterization'],
        ['Assembly Tool:', 'Unicycler'],
    ]
    
    info_table = Table(info_data, colWidths=[2*inch, 4*inch])
    info_table.setStyle(TableStyle([
        ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
        ('FONTNAME', (0, 0), (0, -1), 'Helvetica-Bold'),
        ('FONTSIZE', (0, 0), (-1, -1), 12),
        ('GRID', (0, 0), (-1, -1), 1, colors.black),
        ('BACKGROUND', (0, 0), (0, -1), colors.lightgrey),
    ]))
    story.append(info_table)
    story.append(Spacer(1, 30))
    
    # Assembly Summary
    story.append(Paragraph('Assembly Summary', heading_style))
    
    summary_data = [
        ['Metric', 'Value'],
        ['Plasmids Found', assembly_data.get('plasmid_count', '0')],
        ['Chromosomes Found', assembly_data.get('chromosome_count', '0')],
        ['Total Plasmid Sequences', str(plasmid_count)],
        ['Total Plasmid Length', f'{{plasmid_length:,}} bp' if plasmid_length > 0 else '0 bp'],
    ]
    
    summary_table = Table(summary_data, colWidths=[3*inch, 2*inch])
    summary_table.setStyle(TableStyle([
        ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
        ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
        ('FONTSIZE', (0, 0), (-1, -1), 11),
        ('GRID', (0, 0), (-1, -1), 1, colors.black),
        ('BACKGROUND', (0, 0), (-1, 0), colors.lightblue),
    ]))
    story.append(summary_table)
    story.append(Spacer(1, 20))
    
    # Detailed plasmid information
    if assembly_data.get('plasmids'):
        story.append(Paragraph('Plasmid Details', heading_style))
        
        plasmid_details = [['Contig', 'Size (bp)', 'Classification']]
        for plasmid in assembly_data['plasmids']:
            # Parse plasmid line
            parts = plasmid.split(': ')
            if len(parts) >= 2:
                contig_info = parts[0].strip()
                size_info = parts[1].strip()
                
                # Extract contig number/name
                contig = contig_info.split()[0] if contig_info.split() else 'Unknown'
                
                # Extract size 
                if 'bp' in size_info:
                    size = size_info.split('bp')[0].split()[-1] if size_info.split('bp')[0].split() else 'Unknown'
                else:
                    size = 'Unknown'
                    
                # Extract classification
                if '[' in size_info and ']' in size_info:
                    classification = size_info.split('[')[1].split(']')[0] if '[' in size_info else 'Unknown'
                    # Include copy number if present
                    if '(~' in size_info:
                        copy_info = size_info.split('(~')[1].split(')')[0] if '(~' in size_info else ''
                        classification += f' (~{{copy_info}})'
                else:
                    classification = 'Unknown'
                    
                plasmid_details.append([contig, size, classification])
        
        plasmid_table = Table(plasmid_details, colWidths=[1*inch, 1.5*inch, 3.5*inch])
        plasmid_table.setStyle(TableStyle([
            ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, -1), 10),
            ('GRID', (0, 0), (-1, -1), 1, colors.black),
            ('BACKGROUND', (0, 0), (-1, 0), colors.lightgreen),
        ]))
        story.append(plasmid_table)
        story.append(Spacer(1, 20))
    
    # BLAST Results
    story.append(PageBreak())
    story.append(Paragraph('BLAST Analysis Results', heading_style))
    
    if blast_data:
        blast_table_data = [['Reference Plasmid', 'Identity (%)', 'Alignment Length', 'E-value', 'Bit Score']]
        
        # Sort by bit score (best matches first)
        blast_data_sorted = sorted(blast_data, key=lambda x: x['bitscore'], reverse=True)
        
        for result in blast_data_sorted[:10]:  # Top 10 matches
            blast_table_data.append([
                result['subject'],
                f\"{{result['identity']:.1f}}\",
                str(result['length']),
                f\"{{result['evalue']:.2e}}\",
                f\"{{result['bitscore']:.1f}}\"
            ])
        
        blast_table = Table(blast_table_data, colWidths=[2*inch, 1*inch, 1*inch, 1*inch, 1*inch])
        blast_table.setStyle(TableStyle([
            ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
            ('ALIGN', (1, 1), (-1, -1), 'RIGHT'),  # Right-align numeric columns
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, -1), 9),
            ('GRID', (0, 0), (-1, -1), 1, colors.black),
            ('BACKGROUND', (0, 0), (-1, 0), colors.lightyellow),
        ]))
        story.append(blast_table)
    else:
        story.append(Paragraph('No BLAST results available or no significant matches found.', styles['Normal']))
    
    story.append(Spacer(1, 20))
    
    # MLST Results
    story.append(Paragraph('MLST Typing Results', heading_style))
    mlst_text = mlst_data.get('result', 'No MLST data available')
    story.append(Paragraph(f'MLST Result: {{mlst_text}}', styles['Normal']))
    
    story.append(Spacer(1, 30))
    
    # Footer
    story.append(Paragraph('Report generated by CtGAP 2.0 - Automated bacterial genome analysis pipeline', 
                          styles['Normal']))
    
    # Build the PDF
    doc.build(story)
    print(f'PDF report created: {{output_path}}')

# Main execution
try:
    # Parse input data
    assembly_data = parse_assembly_stats('{input.assembly_stats}')
    blast_data = parse_blast_results('{input.blast_results}')
    mlst_data = parse_mlst_results('{input.mlst_results}')
    plasmid_count, plasmid_length = count_plasmid_sequences('{input.plasmids_fasta}')
    
    # Create output directory
    os.makedirs(os.path.dirname('{output.pdf_report}'), exist_ok=True)
    
    # Generate PDF
    create_pdf_report('{output.pdf_report}', assembly_data, blast_data, mlst_data, plasmid_count, plasmid_length)
    
    print('Plasmid PDF report generated successfully!')
    
except Exception as e:
    print(f'Error generating PDF report: {{e}}', file=sys.stderr)
    import traceback
    traceback.print_exc()
    sys.exit(1)
" > {log} 2>&1
    """

# -----------------------------------------------------------------------------
# 6. BLAST PLASMID
# -----------------------------------------------------------------------------
rule blast_plasmid:
    input:
        plasmid = rules.extract_main_plasmid.output.main_plasmid,
        blastdb = "resources/references/ct/20_Ct_plasmids.nhr"
    output:
        blast = OUTDIR / "{sample}" / "denovo" / "blast" / "blast.plasmid.tab",
        status = OUTDIR / "status" / "denovo.blastn.plasmid.{sample}.txt"
    threads: 12
    conda: "../envs/misc.yaml"
    log: OUTDIR / "{sample}" / "log" / "denovo.blast.plasmid.{sample}.log"
    shell: r"""
        mkdir -p "$(dirname {output.blast})"
        
        # Check if plasmid exists
        if [ -s {input.plasmid} ]; then
            # Create working directory for temp files
            workdir=$(dirname {output.blast})
            tempfile="$workdir/{wildcards.sample}_plasmid_renamed.tmp"
            
            # Rename contig header for sample identification
            awk '
            /^>/ {{
                if (NR==1) {{
                    print ">{wildcards.sample}_plasmid"
                    next
                }} else {{
                    exit
                }}
            }}
            {{
                print $0
            }}
            ' {input.plasmid} > "$tempfile"

            # Run BLAST
            blastn \
              -query "$tempfile" \
              -db resources/references/ct/20_Ct_plasmids \
              -max_target_seqs 5 \
              -outfmt 6 \
              -out {output.blast} \
              > {log} 2>&1
              
            # Clean up
            rm -f "$tempfile"
        else
            echo "No plasmid to BLAST" > {log}
            touch {output.blast}
        fi
        
        touch {output.status}
    """

# -----------------------------------------------------------------------------
# 7. MLST PLASMID
# -----------------------------------------------------------------------------
rule mlst_plasmid:
    input:
        plasmid = rules.extract_main_plasmid.output.main_plasmid
    output:
        mlst = OUTDIR / "{sample}" / "denovo" / "mlst" / "{sample}.genome.plasmid.mlst.txt",
        status = OUTDIR / "status" / "denovo.mlst.plasmid.{sample}.txt"
    conda: "../envs/mlst.yaml"
    log: OUTDIR / "{sample}" / "log" / "denovo.mlst.plasmid.{sample}.log"
    benchmark: OUTDIR / "{sample}" / "benchmark" / "denovo.mlst.plasmid.{sample}.txt"
    shell: r"""
        mkdir -p "$(dirname {output.mlst})"
        
        if [ -s {input.plasmid} ]; then
            echo -e "plasmid\n" > {log}
            claMLST search \
            resources/references/ct/pubmlst/plasmid \
            {input.plasmid} > {output.mlst} 2>> {log}
        else
            echo "No plasmid for MLST typing" > {log}
            touch {output.mlst}
        fi

        touch {output.status}
    """

# -----------------------------------------------------------------------------
# 8. COLLATE BLAST PLASMID
# -----------------------------------------------------------------------------
rule denovo_collate_blast_plasmid:
    input:
        status_blast = expand(OUTDIR / "status" / "denovo.blastn.plasmid.{sample}.txt", sample = SAMPLES),
    output:
        tsv = OUTDIR / "denovo.plasmid.blast.tsv",
        summary = OUTDIR / "denovo.plasmid.blast.summary.tsv",
        status = OUTDIR / "status" / "denovo.plasmid.collate.blast.txt",
    params:
        outdir = OUTDIR,
        pattern = "**/denovo/blast/blast.plasmid.tab",
    threads: 1
    run:
        import pandas as pd
        from pathlib import Path

        # Collect all BLAST results
        all_files = list(Path(params.outdir).glob(params.pattern))

        # Write raw results
        with open(output.tsv, 'w') as out:
            out.write("query\tsubject\tpident\tlength\tmismatch\tgapopen\tquery_start\tquery_end\tsubject_start\tsubject_end\tevalue\tbitscore\n")
            for f in all_files:
                with open(f) as infile:
                    out.write(infile.read())

        # Create summary with cumulative bitscores
        if all_files:
            df = pd.read_csv(output.tsv, sep='\t')
            if not df.empty:
                summary = df.groupby(['query', 'subject']).agg({
                    'bitscore': 'sum',  # Total bitscore
                    'pident': 'mean',   # Average percent identity
                    'length': 'sum'     # Total alignment length
                }).reset_index()
                summary = summary.sort_values(['query', 'bitscore'], ascending=[True, False])
                summary.to_csv(output.summary, sep='\t', index=False)
            else:
                with open(output.summary, 'w') as out:
                    out.write("query\tsubject\tbitscore\tpident\tlength\n")
        else:
            with open(output.summary, 'w') as out:
                out.write("query\tsubject\tbitscore\tpident\tlength\n")

        Path(output.status).touch()

# -----------------------------------------------------------------------------
# 9. COLLATE MLST PLASMID
# -----------------------------------------------------------------------------
rule denovo_collate_mlst_plasmid:
    input:
        plasmid = expand(OUTDIR / "{sample}" / "denovo" / "mlst" / "{sample}.genome.plasmid.mlst.txt", sample = SAMPLES),
    output:
        plasmid = OUTDIR / "denovo.mlst.plasmid.results.tsv",
        status = OUTDIR / "status" / "denovo.mlst.plasmid.collate.txt"
    conda: "../envs/misc.yaml"
    threads: 1
    shell: r"""
        csvtk concat {input.plasmid} -o {output.plasmid}
        touch {output.status}
    """

# -----------------------------------------------------------------------------
# 10. UNICYCLER SUMMARY REPORT
# -----------------------------------------------------------------------------
rule unicycler_summary:
    input:
        expand(OUTDIR / "{sample}" / "denovo" / "unicycler" / "assembly.fasta", sample=SAMPLES)
    output:
        report = OUTDIR / "reports" / "unicycler_summary.txt"
    shell: r"""
      mkdir -p "$(dirname {output.report})"

      echo "Unicycler Assembly Summary Report" > {output.report}
      echo "Generated: $(date)" >> {output.report}
      echo "=================================" >> {output.report}
      echo "" >> {output.report}

      # Count successes and failures
      success=$(find {OUTDIR}/*/denovo/unicycler -name "UNICYCLER_SUCCESS.txt" 2>/dev/null | wc -l)
      failed=$(find {OUTDIR}/*/denovo/unicycler -name "UNICYCLER_FAILED.txt" 2>/dev/null | wc -l)

      echo "Total samples processed: $(($success + $failed))" >> {output.report}
      echo "Successful assemblies: $success" >> {output.report}
      echo "Failed assemblies: $failed" >> {output.report}
      echo "" >> {output.report}

      if [ $failed -gt 0 ]; then
        echo "Samples with failed assemblies:" >> {output.report}
        find {OUTDIR}/*/denovo/unicycler -name "UNICYCLER_FAILED.txt" 2>/dev/null | \
          sed 's|{OUTDIR}/||' | sed 's|/denovo/unicycler/UNICYCLER_FAILED.txt||' | \
          sed 's/^/  - /' >> {output.report}
        echo "" >> {output.report}
      fi

      if [ $success -gt 0 ]; then
        echo "Successful assemblies with plasmid analysis:" >> {output.report}
        echo "" >> {output.report}
        
        for sample_dir in $(find {OUTDIR}/*/denovo/unicycler -name "UNICYCLER_SUCCESS.txt" 2>/dev/null | sed 's|/UNICYCLER_SUCCESS.txt||'); do
          sample=$(echo $sample_dir | sed 's|{OUTDIR}/||' | sed 's|/denovo/unicycler||')
          stats_file="$sample_dir/assembly_stats.txt"
          plasmid_file="$sample_dir/plasmids_only.fasta"
          
          echo "  Sample: $sample" >> {output.report}
          
          if [ -f "$stats_file" ]; then
            # Extract key info from stats
            plasmid_count=$(grep "Plasmids found:" "$stats_file" | cut -d: -f2 | tr -d ' ' 2>/dev/null || echo "0")
            echo "    Plasmids found: $plasmid_count" >> {output.report}
            
            if [ -f "$plasmid_file" ] && [ -s "$plasmid_file" ]; then
              echo "    Plasmid sizes:" >> {output.report}
              awk '/^>/ {{
                if(seq && header) printf "      %s: %d bp\\n", header, length(seq)
                header = $0; gsub(/^>/, "", header); seq = ""; next
              }} 
              {{seq = seq $0}} 
              END {{if(seq && header) printf "      %s: %d bp\\n", header, length(seq)}}' "$plasmid_file" >> {output.report}
            else
              echo "    No plasmids extracted" >> {output.report}
            fi
          fi
          echo "" >> {output.report}
        done
      fi

      echo "Check individual files:" >> {output.report}
      echo "  Assembly: {OUTDIR}/<sample>/denovo/unicycler/assembly.fasta" >> {output.report}
      echo "  Plasmids: {OUTDIR}/<sample>/denovo/unicycler/plasmids_only.fasta" >> {output.report}
      echo "  Stats: {OUTDIR}/<sample>/denovo/unicycler/assembly_stats.txt" >> {output.report}
      echo "  Logs: {OUTDIR}/<sample>/log/denovo.unicycler.<sample>.log" >> {output.report}
      echo "  PDF Reports: {OUTDIR}/<sample>/reports/<sample>_plasmid_report.pdf" >> {output.report}
    """

# -----------------------------------------------------------------------------
# 11. COMPREHENSIVE PLASMID ANALYSIS SUMMARY
# -----------------------------------------------------------------------------
rule plasmid_analysis_summary:
    input:
        blast_summary = rules.denovo_collate_blast_plasmid.output.summary,
        mlst_results = rules.denovo_collate_mlst_plasmid.output.plasmid,
        unicycler_report = rules.unicycler_summary.output.report,
        pdf_reports = expand(OUTDIR / "{sample}" / "reports" / "{sample}_plasmid_report.pdf", sample=SAMPLES)
    output:
        comprehensive_report = OUTDIR / "reports" / "comprehensive_plasmid_analysis.txt"
    shell: r"""
        mkdir -p "$(dirname {output.comprehensive_report})"

        echo "COMPREHENSIVE PLASMID ANALYSIS REPORT" > {output.comprehensive_report}
        echo "====================================" >> {output.comprehensive_report}
        echo "Generated: $(date)" >> {output.comprehensive_report}
        echo "" >> {output.comprehensive_report}

        echo "This report summarizes all plasmid-related analyses:" >> {output.comprehensive_report}
        echo "1. Plasmid assembly (Unicycler)" >> {output.comprehensive_report}
        echo "2. Plasmid/chromosome separation" >> {output.comprehensive_report}
        echo "3. Main plasmid extraction" >> {output.comprehensive_report}
        echo "4. Plasmid BLAST detection (vs 20 reference plasmids)" >> {output.comprehensive_report}
        echo "5. Plasmid MLST typing" >> {output.comprehensive_report}
        echo "6. Individual PDF reports per sample" >> {output.comprehensive_report}
        echo "" >> {output.comprehensive_report}

        echo "=== INDIVIDUAL PDF REPORTS ===" >> {output.comprehensive_report}
        echo "" >> {output.comprehensive_report}
        echo "Sample-specific PDF reports generated:" >> {output.comprehensive_report}
        for pdf_file in {input.pdf_reports}; do
            sample=$(basename "$pdf_file" | sed 's/_plasmid_report.pdf//')
            echo "  $sample: $pdf_file" >> {output.comprehensive_report}
        done
        echo "" >> {output.comprehensive_report}

        echo "=== DETAILED REPORTS ===" >> {output.comprehensive_report}
        echo "" >> {output.comprehensive_report}
        echo "For detailed results, see:" >> {output.comprehensive_report}
        echo "• BLAST results: {OUTDIR}/denovo.plasmid.blast.summary.tsv" >> {output.comprehensive_report}
        echo "• MLST results: {OUTDIR}/denovo.mlst.plasmid.results.tsv" >> {output.comprehensive_report}
        echo "• Unicycler report: {input.unicycler_report}" >> {output.comprehensive_report}
        echo "" >> {output.comprehensive_report}

        echo "=== UNICYCLER ASSEMBLY SUMMARY ===" >> {output.comprehensive_report}
        tail -n +4 {input.unicycler_report} >> {output.comprehensive_report}
    """
