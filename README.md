# CtGAP 2.0

***Chlamydia trachomatis* Genome Assembly Pipeline**

CtGAP assembles *C. trachomatis* genomes from Illumina paired-end reads, performs strain typing, OmpA genotyping,MLST, detects plasmids (and types them), and builds phylogenetic trees.

---

## What's Included

CtGAP has all databases and references included:

- Human reference genome (GRCh38) for host read removal
- Kraken2 database for taxonomic classification
- *C. trachomatis* reference genomes
- OmpA database (26 genotypes, and growing)
- MLST schemes (Chlamydiales and *C. trachomatis*)
- 20 plasmid reference sequences

*Large databases (GRCh38, Kraken2) are downloaded automatically during setup.

**No additional downloads required.**

---

### Dependencies

- Conda ([install guide](#installing-conda))
- Snakemake ≥7.0 ([install guide](#installing-snakemake))

---

## Installation

```bash
git clone https://github.com/D-Dean-Lab/CtGAP
cd ctgap
ctgap setup
```

setup will:
1. Download required databases (GRCh38 + Kraken2)
2. Create Conda environments

---

## Quick Start

### 1. Add reads

Copy paired-end FASTQ files to `input/`:

**File naming must follow this pattern:** `{sample}_R1.fastq.gz` and `{sample}_R2.fastq.gz`

| Valid ✓ | Invalid ✗ |
|---------|-----------|
| `CT001_R1.fastq.gz` | `CT001_1.fastq.gz` |
| `patient_A_R1.fastq.gz` | `patient_A_.R1.fastq.gz` |
| `2024_sample01_R1.fastq.gz` | `2024_sample01_R1.fq.gz` |

### 2. Run

```bash
ctgap run
```

### 3. Get results (truncated list; you get way more)

| File | Description |
|------|-------------|
| `output/{sample}/best/assembly.fasta` | Final genome assembly |
| `output/ct.tree` | Phylogenetic tree (Newick) |
| `output/best.ompA_genovar.blast.tsv` | OmpA genovar typing |
| `output/best.mlst.ct.results.tsv` | MLST sequence types |
| `output/{sample}/reports/{sample}_plasmid_report.pdf` | Plasmid report |

---

## Assembly Modes

```bash
ctgap run                        # Auto (default)
ctgap run --denovo               # De novo only
ctgap run --ref-guided [REF]     # Reference-guided with choice reference
```

| Mode | When to use |
|------|-------------|
| **Auto** | Most cases — runs both methods, picks best per sample |
| **De novo** | Divergent strains, avoiding reference bias |
| **Reference-guided** | Better contiguity with known reference |

Available references: any genome in `resources/references/ct/individual/`, or `plurality` for consensus refseq.

### Run time estimates

| Samples | Mode | Time* |
|---------|------|-------|
| 1 | Auto | 30–60 min |
| 10 | Auto | 2–4 hours |
| 50 | Auto | 8–12 hours |

*based on a 48-core system, 64 GB RAM. Varies with read depth.

---

## Optional Configuration

Edit `config/config.yaml` to adjust thread counts for your system:

---

## Output (truncated list; you get way more)

```
output/
├── {sample}/
│   ├── best/
│   │   └── assembly.fasta              # Final assembly
│   ├── denovo/unicycler/
│   │   └── plasmids_only.fasta         # Detected plasmids
│   └── reports/
│       └── {sample}_plasmid_report.pdf # Plasmid report
│
├── ct.tree                             # Phylogenetic tree
├── best.ompA_genovar.blast.tsv         # OmpA results (all samples)
├── best.mlst.ct.results.tsv            # MLST results (all samples)
├── denovo.mlst.plasmid.results.tsv     # Plasmid MLST (all samples)
│
└── reports/
    ├── assembly_methods_summary.txt    # Method selected per sample
    └── comprehensive_plasmid_analysis.txt
```

---

## Troubleshooting

### "No samples were detected"

Files aren't named correctly. Check the pattern:

```bash
ls input/
# Must be {sample}_R1.fastq.gz and {sample}_R2.fastq.gz
```

### Pipeline failed mid-run?

Just re-run. It resumes automatically:

```bash
ctgap run
```

### View logs

```bash
cat output/{sample}/log/*.log   # Per-sample logs
cat .snakemake/log/*.log        # Snakemake logs
```

---

## Feature Summary

| Feature | Description |
|---------|-------------|
| **Automated assembly selection** | Runs de novo and reference-guided, scores each (N50, contiguity, GC%), picks best |
| **Host read removal** | Depletes human reads using HoCort + Bowtie2 before assembly |
| **Taxonomic filtering** | Extracts Chlamydiales reads via Kraken2 for cleaner assemblies |
| **Hybrid assembly** | Shovill/SPAdes assembly with RagTag scaffolding and Gap2Seq gap-filling |
| **OmpA genovar typing** | BLAST against 20 reference OmpA sequences with secondary typing |
| **MLST typing** | Dual-scheme typing (Chlamydiales + *C. trachomatis*) via pyMLST |
| **Plasmid detection** | Unicycler-based assembly, chromosome/plasmid separation, copy number estimation |
| **Plasmid typing** | BLAST against 20 plasmid references + plasmid MLST |
| **Phylogenetics** | SKA2 core genome alignment + IQ-TREE with bootstrap support |
| **PDF reports** | Per-sample plasmid analysis reports |
| **QC integration** | QUAST assembly stats, coverage analysis, MultiQC summary |
| **Batteries-included** | All databases pre-configured — GRCh38, Kraken2, MLST, BLAST, references |
| **Resumable** | Automatic resume from failure point |
| **Scalable** | Processes 1 to 100+ samples with configurable parallelism |

---

### Installing Conda

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc
```

### Installing Snakemake

```bash
conda install -n base -c conda-forge mamba
mamba create -n snakemake -c conda-forge -c bioconda snakemake
conda activate snakemake
```

---
