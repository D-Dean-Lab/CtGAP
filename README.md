# CtGAP 2.0

***Chlamydia trachomatis* Genome Assembly Pipeline**

CtGAP assembles *C. trachomatis* genomes from Illumina paired-end reads, performs comprehensive typing (ompA genotype, MLST, whole-genome strain typing), detects and types plasmids, annotates genomes, detects mixed infections, and builds phylogenetic trees.

---

## What's New in v2.0 (Highlight)

- **Automatic assembly selection** — Runs both de novo and reference-guided, picks best per sample
- **Whole-genome strain typing** — Mash + FastANI against 21 reference genomes
- **Mixed strain detection** — Flags potential co-infections via heterozygous site analysis
- **Genome annotation** — Bakta annotation (optional)
- **Dedicated plasmid module** — Reference-guided extraction + Unicycler assembly
- **Comprehensive PDF reports** — Per-sample reports with all results in one document
- **Assembly filtering** — Removes non-Ct contigs automatically

---

## What's Included

CtGAP comes with all databases and references:

- Human reference genome (GRCh38) for host read removal
- Kraken2 database for taxonomic classification
- *C. trachomatis* reference genomes
- OmpA database (26 genotypes, and growing)
- MLST schemes (Chlamydiales and *C. trachomatis*)
- 20 plasmid reference sequences
- Bakta database for annotation

*Large databases (GRCh38, Kraken2, Bakta) are downloaded automatically during setup.*

**No additional downloads required.**

---

## Dependencies

- Conda ([install guide](#installing-conda))
- Snakemake ≥7.0 ([install guide](#installing-snakemake))

---

## Installation

```bash
git clone https://github.com/D-Dean-Lab/CtGAP
cd CtGAP
./ctgap setup
```

Setup will:
1. Download required databases (GRCh38, Kraken2, Bakta)
2. Create Conda environments
3. Build BLAST databases

---

## Quick Start

### 1. Add reads

Copy paired-end FASTQ files to `input/`:

**File naming must follow this pattern:** `{sample}_R1.fastq.gz` and `{sample}_R2.fastq.gz`

| Valid | Invalid |
|-------|---------|
| `CT001_R1.fastq.gz` | `CT001_1.fastq.gz` |
| `patient_A_R1.fastq.gz` | `patient_A_.R1.fastq.gz` |
| `2024_sample01_R1.fastq.gz` | `2024_sample01_R1.fq.gz` |

### 2. Run

```bash
./ctgap run
```

### 3. Get results (truncated list; you get way more)

After completion, run `./ctgap clean` to organize outputs:

```
output/
├── 01_per_sample/{sample}/
│   ├── best/assembly.fasta           # Final genome assembly
│   ├── plasmid/{sample}.plasmid.fasta
│   └── reports/{sample}_ctgap_report.pdf
│
├── 02_reports/
│   ├── ctgap_all_samples_summary.csv  # All samples, one row each
│   ├── ctgap_batch_summary.txt
│   └── multiqc_report.html
│
├── 03_results/
│   ├── typing/best.ompA_genovar.blast.tsv
│   ├── typing/best.mlst.ct.results.tsv
│   └── plasmid/plasmid_summary.tsv
│
└── 04_phylogeny/ct.tree               # Newick format
```

---

## Assembly Modes

```bash
./ctgap run                        # Auto (default) — recommended
./ctgap run --denovo               # De novo only
./ctgap run --ref-guided           # Reference-guided (plurality consensus)
./ctgap run --ref-guided E_bour    # Reference-guided with specific strain
```

| Mode | When to use |
|------|-------------|
| **Auto** | Most cases — runs both methods, scores each, picks best per sample |
| **De novo** | Divergent strains, avoiding reference bias |
| **Reference-guided** | Better contiguity when reference is known |

**Available references:** `plurality` (consensus), or any genome in `resources/references/ct/individual/`

---

## Assembly Selection Scoring

When using **auto** mode, CtGAP scores each assembly:

| Metric | Weight | Criteria |
|--------|--------|----------|
| N50 | 40% | ≥95% of expected genome = full score |
| Contig count | 30% | ≤2 contigs = full score |
| Length accuracy | 20% | Deviation from 1.04 Mb expected |
| GC content | 10% | Optimal: 40-42% |

The assembly with the higher score is selected. Ties favor de novo (less reference bias).

---

## Typing Results

| Analysis | Method | Output |
|----------|--------|--------|
| **ompA Genotype** | BLAST (+ secondary disambiguation) | `best.ompA_genovar.blast.tsv` |
| **MLST** | Chlamydiales + *C. trachomatis* schemes | `best.mlst.*.results.tsv` |
| **Ct Strain** | Mash distance + FastANI | `ct_typing_all_samples.csv` |
| **Plasmid** | Dedicated assembly + BLAST + MLST | `plasmid_summary.tsv` |

---

## Optional Features (see config for editable parameters)

Enable/disable in `config/config.yaml`:

```yaml
# Genome annotation (Bakta)
annotation:
  enabled: true        # Set false to skip

# Ct strain typing (Mash + FastANI)
ct_typing:
  enabled: true        # Set false to skip

# Mixed strain detection
mixed_detection:
  enabled: true        # Set false to skip
  het_threshold: 50    # Sites to flag as mixed

# Assembly filtering (remove non-Ct contigs)
assembly_filter:
  enabled: true
  min_coverage: 30     # % of contig aligned to Ct reference
```

---

## PDF Reports

Each sample gets a comprehensive PDF report containing:

- **Key Findings** — Genotype, strain, plasmid status at a glance
- **Assembly Metrics** — Length, N50, contigs, gaps in main contig, GC%
- **Typing Results** — ompA, MLST, Ct strain assignment
- **Plasmid Analysis** — Detection, size, circularity, BLAST hits
- **QC Metrics** — Read stats, coverage, Kraken2 classification
- **Mixed Strain Status** — Heterozygous sites, interpretation

---

## Troubleshooting

### "No samples were detected"

Files aren't named correctly:
```bash
ls input/
# Must be {sample}_R1.fastq.gz and {sample}_R2.fastq.gz
```

### Pipeline failed mid-run?

Just re-run — it resumes automatically:
```bash
./ctgap run
```

### Check logs

```bash
cat output/{sample}/log/*.log   # Per-sample logs
cat .snakemake/log/*.log        # Snakemake logs
```

### Low assembly quality?

Check coverage in the PDF report. Low input depth (<10x) affects assembly contiguity.

---

## Run Time Estimates

| Samples | Mode | Time* |
|---------|------|-------|
| 1 | Auto | 30–60 min |
| 10 | Auto | 2–4 hours |
| 50 | Auto | 8–12 hours |

*Based on 32-core system, 64 GB RAM. Varies with read depth and optional modules enabled.

---

## Feature Summary

| Feature | Description |
|---------|-------------|
| **Auto assembly selection** | Scores de novo and ref-guided, picks best per sample |
| **Host read removal** | HoCoRT + Bowtie2 depletion of human reads |
| **Taxonomic filtering** | Kraken2 extraction of Chlamydiales reads |
| **Assembly** | Shovill/SPAdes + RagTag scaffolding + Gap2Seq |
| **Assembly filtering** | Removes non-Ct contigs via minimap2 |
| **ompA genotyping** | BLAST with secondary disambiguation (A/B/Ba/C) |
| **MLST** | Dual-scheme (Chlamydiales + *C. trachomatis*) |
| **Ct strain typing** | Whole-genome Mash + FastANI against 21 references |
| **Mixed detection** | Heterozygous site analysis for co-infections |
| **Plasmid module** | Reference-guided extraction + Unicycler + BLAST + MLST |
| **Annotation** | Bakta genome annotation (optional) |
| **Phylogenetics** | SKA2 alignment + IQ-TREE |
| **PDF reports** | Comprehensive per-sample reports |
| **QC integration** | QUAST, MultiQC, coverage analysis |
| **Resumable** | Automatic checkpoint recovery |

---

## Installing Conda

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc
```

## Installing Snakemake

```bash
conda install -n base -c conda-forge mamba
mamba create -n snakemake -c conda-forge -c bioconda snakemake
conda activate snakemake
```

---

## Dependencies

[fastp](https://github.com/OpenGene/fastp) ·
[HoCoRT](https://github.com/ignasrum/hocort) ·
[Bowtie2](https://github.com/BenLangmead/bowtie2) ·
[SAMtools](https://github.com/samtools/samtools) ·
[BCFtools](https://github.com/samtools/bcftools) ·
[Kraken2](https://github.com/DerrickWood/kraken2) ·
[KrakenTools](https://github.com/jenniferlu717/KrakenTools) ·
[Shovill](https://github.com/tseemann/shovill) ·
[SPAdes](https://github.com/ablab/spades) ·
[Unicycler](https://github.com/rrwick/Unicycler) ·
[RagTag](https://github.com/malonge/RagTag) ·
[Gap2Seq](https://github.com/rikuturkki/Gap2Seq) ·
[Pilon](https://github.com/broadinstitute/pilon) ·
[minimap2](https://github.com/lh3/minimap2) ·
[BLAST](https://blast.ncbi.nlm.nih.gov) ·
[claMLST](https://github.com/bvalot/pyMLST) ·
[Mash](https://github.com/marbl/Mash) ·
[FastANI](https://github.com/ParBLiSS/FastANI) ·
[Bakta](https://github.com/oschwengers/bakta) ·
[QUAST](https://github.com/ablab/quast) ·
[SKA2](https://github.com/bacpop/ska.rust) ·
[IQ-TREE](https://github.com/iqtree/iqtree2) ·
[MultiQC](https://github.com/MultiQC/MultiQC) ·
[ReportLab](https://www.reportlab.com)

---

## Citation

If you use CtGAP in your research, please cite


---

## License

MIT License - see [LICENSE](LICENSE) for details.
