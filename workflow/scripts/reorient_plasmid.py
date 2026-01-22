#!/usr/bin/env python3
"""
Re-orient circular plasmid to match reference orientation and start position.

Circular plasmids may be assembled in reverse complement or with a different
start position. This causes BLAST to report split alignments.
"""

import subprocess
import sys
from pathlib import Path


def reverse_complement(seq):
    """Return reverse complement of a DNA sequence."""
    complement = dict(A='T', T='A', G='C', C='G', a='t', t='a', g='c', c='g', N='N', n='n')
    return ''.join(complement.get(base, base) for base in reversed(seq))


def rotate_sequence(seq, position):
    """Rotate sequence so it starts at given position."""
    if position <= 0 or position >= len(seq):
        return seq
    return seq[position:] + seq[:position]


def parse_fasta(filename):
    """Parse FASTA file, return list of (header, sequence) tuples."""
    sequences = []
    with open(filename) as f:
        header = ""
        seq = ""
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header and seq:
                    sequences.append((header, seq))
                header = line[1:]
                seq = ""
            else:
                seq += line
        if header and seq:
            sequences.append((header, seq))
    return sequences


def write_fasta(filename, header, seq):
    """Write sequence to FASTA file with 80-char line wrapping."""
    with open(filename, 'w') as f:
        f.write(">" + header + "\n")
        for i in range(0, len(seq), 80):
            f.write(seq[i:i+80] + "\n")


def main():
    # Get inputs from snakemake
    plasmid_file = snakemake.input.plasmid
    reference_file = snakemake.input.reference
    output_file = snakemake.output.reoriented
    log_file = snakemake.output.orientation_log

    # Parse input plasmid
    plasmid_seqs = parse_fasta(plasmid_file)

    if not plasmid_seqs or "no_plasmid" in plasmid_seqs[0][0] or "no_valid" in plasmid_seqs[0][0]:
        # No valid plasmid, just copy
        with open(output_file, 'w') as f:
            for header, seq in plasmid_seqs:
                f.write(">" + header + "\n")
                for i in range(0, len(seq), 80):
                    f.write(seq[i:i+80] + "\n")
        with open(log_file, 'w') as f:
            f.write("No valid plasmid to reorient\n")
        return

    header, seq = plasmid_seqs[0]
    original_len = len(seq)

    # Run minimap2 to align plasmid to reference
    try:
        result = subprocess.run(
            ["minimap2", "-cx", "asm5", reference_file, plasmid_file],
            capture_output=True, text=True, timeout=60
        )
        paf_lines = result.stdout.strip().split("\n")
    except Exception as e:
        # If minimap2 fails, just copy original
        write_fasta(output_file, header, seq)
        with open(log_file, 'w') as f:
            f.write("minimap2 failed: " + str(e) + "\nKept original orientation\n")
        return

    # Parse PAF to find best alignment
    # PAF format: qname qlen qstart qend strand tname tlen tstart tend matches alnlen mapq
    best_alignment = None
    best_matches = 0

    for line in paf_lines:
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) >= 12:
            matches = int(fields[9])
            if matches > best_matches:
                best_matches = matches
                best_alignment = {
                    'qstart': int(fields[2]),
                    'qend': int(fields[3]),
                    'strand': fields[4],
                    'tstart': int(fields[7]),
                    'tend': int(fields[8]),
                    'matches': matches,
                    'alnlen': int(fields[10])
                }

    if not best_alignment:
        # No alignment found, keep original
        write_fasta(output_file, header, seq)
        with open(log_file, 'w') as f:
            f.write("No alignment found\nKept original orientation\n")
        return

    # Determine if reverse complement is needed
    needs_revcomp = best_alignment['strand'] == '-'

    # Apply reverse complement if needed
    if needs_revcomp:
        seq = reverse_complement(seq)
        # Estimate rotation based on original alignment
        rotation_pos = original_len - best_alignment['qend']
    else:
        rotation_pos = best_alignment['qstart'] - best_alignment['tstart']

    # Normalize rotation position
    if rotation_pos < 0:
        rotation_pos = original_len + rotation_pos
    rotation_pos = rotation_pos % original_len

    # Rotate sequence
    if rotation_pos > 0:
        seq = rotate_sequence(seq, rotation_pos)

    # Write output
    write_fasta(output_file, header + " reoriented", seq)

    # Write orientation log
    with open(log_file, 'w') as f:
        f.write("Original length: " + str(original_len) + " bp\n")
        f.write("Alignment matches: " + str(best_alignment['matches']) + " / " + str(best_alignment['alnlen']) + "\n")
        f.write("Strand: " + best_alignment['strand'] + "\n")
        f.write("Reverse complemented: " + ("Yes" if needs_revcomp else "No") + "\n")
        f.write("Rotation applied: " + str(rotation_pos) + " bp\n")
        f.write("Final length: " + str(len(seq)) + " bp\n")

    print("Plasmid reoriented: revcomp=" + str(needs_revcomp) + ", rotation=" + str(rotation_pos))


if __name__ == "__main__":
    main()
