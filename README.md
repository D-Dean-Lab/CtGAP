# CtGAP - Chlamydia trachomatis Genome Assembly Pipeline

### Install

1. `git clone` this repo:

```
git clone https://github.com/ammaraziz/ctgap
```

2. Install `miniconda` or preferrably `mamba`
	- https://mamba.readthedocs.io/en/latest/installation.html
	- https://docs.conda.io/projects/miniconda/en/latest/

3. Install `snakemake`:
```
mamba install -c bioconda snakemake 
```

Manually install rust/scrubby. Make sure it is the newest release (v0.3.0 updated 07/01/2024)
```
git clone https://github.com/esteinig/scrubby
cd scrubby && cargo build --release
./target/release/scrubby --help
```
4. Install `kraken2`:
```
./install_kraken2.sh $KRAKEN2_DIR
cp $KRAKEN2_DIR/kraken2{,-build,-inspect} $HOME/bin
```
(Replace `$KRAKEN2_DIR` above with the directory where you want to install Kraken 2's programs/scripts.)

5. Download the human genome, rename to `resources/grch38.fasta`
	- [From NCBI](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/)

6. Download one of the kraken dbs with bacterial genomes, rename to `resources/standardDB`:
	- https://benlangmead.github.io/aws-indexes/k2

7. Done - The pipeline will handle the dependencies internally. Ensure all 
downloaded packages and your kraken2 directory are on your path. 

### Usage

1. Create a folder `ctgap/input/`
2. Add your fastq.gz files in `ctgap/input/`. 

	- Ensure they're named as follows: `{sample_name}_{direction}.fastq.gz`. 
        - eg `SRR12345_R1.fastq.gz` and `SRR12345_R2.fastq.gz`.

3. In `ctgap/` folder run the pipeline:
```
snakemake -j 8 --use-conda -k
```

- `-j 8` specifies the number of threads to use in total.
- `--use-conda` tells snakemake to install the dependencies.
- `-k` tells snakemake to keep going if a sample fails.

### Dependencies

- Snakemake
- Spades
- Shovill
- Bowtie2
- Samtools
- fastp
- bbmap (bbnorm)
- kraken2
- multiqc
- blast+

### Output

TBA

### Cite

Pipeline is created by Shola Olagoke with assistance from Ammar Aziz.
