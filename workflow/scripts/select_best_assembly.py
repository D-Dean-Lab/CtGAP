#!/usr/bin/env python3
"""
Select the best assembly (denovo vs reference-guided) based on quality metrics
"""
import sys
import shutil
from pathlib import Path

def parse_quast_report(quast_file):
    """Parse QUAST report to extract key metrics"""
    stats = {}
    with open(quast_file, 'r') as f:
        for line in f:
            if '\t' in line:
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    key, value = parts[0], parts[1]
                    if key == 'Total length':
                        stats['total_length'] = int(value.replace(',', ''))
                    elif key == 'N50':
                        stats['N50'] = int(value.replace(',', ''))
                    elif key == '# contigs':
                        stats['num_contigs'] = int(value)
                    elif key == 'GC (%)':
                        stats['gc_percent'] = float(value)
    return stats

def calculate_assembly_score(stats, expected_length=1040000):
    """Calculate quality score (0-100 scale)"""
    score = 0
    
    # N50 Score (40 points)
    n50_ratio = stats['N50'] / expected_length
    if n50_ratio >= 0.95:
        score += 40
    elif n50_ratio >= 0.5:
        score += 40 * n50_ratio
    else:
        score += 40 * (n50_ratio ** 2)
    
    # Contig Count Score (30 points)
    num_contigs = stats['num_contigs']
    if num_contigs <= 2:
        score += 30
    elif num_contigs <= 5:
        score += 25
    elif num_contigs <= 10:
        score += 20
    elif num_contigs <= 20:
        score += 10
    else:
        score += max(0, 10 - (num_contigs - 20))
    
    # Length Accuracy Score (20 points)
    length_diff = abs(stats['total_length'] - expected_length)
    if length_diff <= 10000:
        score += 20
    elif length_diff <= 50000:
        score += 20 * (1 - length_diff/50000)
    else:
        score += max(0, 10 - (length_diff - 50000)/10000)
    
    # GC Content Check (10 points)
    gc_content = stats['gc_percent']
    if 40.0 <= gc_content <= 42.0:
        score += 10
    elif 39.0 <= gc_content <= 43.0:
        score += 8
    else:
        score += max(0, 8 - abs(gc_content - 41.0))
    
    return round(score, 2)

def main():
    # Parse input from Snakemake
    denovo_stats_file = snakemake.input.denovo_quast
    refguided_stats_file = snakemake.input.refguided_quast
    denovo_asm = snakemake.input.denovo
    refguided_asm = snakemake.input.refguided
    output_asm = snakemake.output.best
    report_file = snakemake.output.report
    metadata_file = snakemake.output.metadata
    
    # Parse QUAST reports
    denovo_stats = parse_quast_report(denovo_stats_file)
    refguided_stats = parse_quast_report(refguided_stats_file)
    
    # Calculate scores
    denovo_score = calculate_assembly_score(denovo_stats)
    refguided_score = calculate_assembly_score(refguided_stats)
    
    # Select best assembly
    if denovo_score >= refguided_score:
        selected = "denovo"
        selected_asm = denovo_asm
        winner_score = denovo_score
        loser_score = refguided_score
    else:
        selected = "reference-denovo"
        selected_asm = refguided_asm
        winner_score = refguided_score
        loser_score = denovo_score
    
    # Copy best assembly to output
    shutil.copy(selected_asm, output_asm)
    
    # Write method metadata
    with open(metadata_file, 'w') as f:
        f.write(selected)
    
    # Write detailed selection report
    with open(report_file, 'w') as f:
        f.write(f"Assembly Selection Report\n")
        f.write(f"========================\n\n")
        f.write(f"SELECTED: {selected} (Score: {winner_score:.2f})\n\n")
        
        f.write(f"De Novo Assembly:\n")
        f.write(f"  Score: {denovo_score:.2f}\n")
        f.write(f"  N50: {denovo_stats['N50']:,} bp\n")
        f.write(f"  Contigs: {denovo_stats['num_contigs']}\n")
        f.write(f"  Length: {denovo_stats['total_length']:,} bp\n")
        f.write(f"  GC%: {denovo_stats['gc_percent']:.2f}\n\n")
        
        f.write(f"Reference-Guided Assembly:\n")
        f.write(f"  Score: {refguided_score:.2f}\n")
        f.write(f"  N50: {refguided_stats['N50']:,} bp\n")
        f.write(f"  Contigs: {refguided_stats['num_contigs']}\n")
        f.write(f"  Length: {refguided_stats['total_length']:,} bp\n")
        f.write(f"  GC%: {refguided_stats['gc_percent']:.2f}\n\n")
        
        f.write(f"Score Difference: {abs(winner_score - loser_score):.2f} points\n")

if __name__ == "__main__":
    main()
