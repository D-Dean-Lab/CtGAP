# CtGAP - Chlamydia trachomatis Genome Assembly Pipeline

### Install

1. `git clone` this repo:

```
git clone https://github.com/D-Dean-Lab/CtGAP
```

2. Install `miniconda` or preferrably `mamba`
	- https://mamba.readthedocs.io/en/latest/installation.html
	- https://docs.conda.io/projects/miniconda/en/latest/

3. Install `snakemake`:
```
mamba install -c bioconda snakemake 
```

4a. Install rust and cmake if you don't already have them
```
mamba install -c conda-forge rust
mamba install -c conda-forge cmake
```
- Ensure cargo is on your path: `export PATH=$PATH:/path/to/.cargo/bin`

4b. Manually install scrubby. Use v0.3.0 (updated 07/01/2024) from the "empty" branch.
```
git clone https://github.com/esteinig/scrubby --branch empty
cd scrubby && cargo build --release
./target/release/scrubby --help
```
- Ensure scrubby is on your path: `export PATH=$PATH:/path/to/scrubby/target/release`
5. Install `kraken2`:
```
./install_kraken2.sh $KRAKEN2_DIR
cp $KRAKEN2_DIR/kraken2{,-build,-inspect} $HOME/bin
```
- (Replace `$KRAKEN2_DIR` above with the directory where you want to install Kraken 2's programs/scripts.)
- Ensure the kraken2 directory is on your path: `export PATH=$PATH:/path/to/kraken2_dir`

6. Download a host genome, and ensure it gets placed into the `resources/` folder
	- We recommend this human  genome[from NCBI](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/)

7. Download one of the kraken dbs with archaea, eukaryotic, and bacterial genomes, rename to `resources/standardDB`:
	- https://benlangmead.github.io/aws-indexes/k2
		- We suggest: https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_08gb_20240904.tar.gz 

8. Done - The pipeline will handle the dependencies internally. Ensure all 
downloaded packages and your kraken2 directory are on your path. 

### Usage

1. Create a folder `ctgap/input/`
2. Add your fastq.gz files to analyze in `ctgap/input/`. 

	- Ensure they're named as follows: `{sample_name}_{direction}.fastq.gz`. 
        - eg `SRR12345_R1.fastq.gz` and `SRR12345_R2.fastq.gz`.

3. Update ctgap/config/config.yaml as needed for your run. Below are values you can update.
- `indir` is the path to your input folder
- `outdir` is the path to your desired output folder. 
- `mode` will tell the pipeline what type of assembly to perform 
- `reference` will be the reference sequence your samples will be compared to. If you are running denovo assembly or don't want to use a reference, set this to "reference_24"
- `hostReference` will be the genome of the host organism you downloaded above in step 6. Add its full filename after `resource/`
4. In `ctgap/` folder run the pipeline:
```
snakemake -j 8 --use-conda -k
```

- `-j 8` specifies the number of threads to use in total. You can change this number based on your needs.
- `--use-conda` tells snakemake to install the dependencies.
- `-k` tells snakemake to keep going if a sample fails.

### Dependencies

- Snakemake
- Spades
- Shovill
- Bowtie2
- Samtools=1.19
- Bcftools
- Bedtools
- Mummer
- fastp
- scrubby
- kraken2
- minimap2
- multiqc
- quast
- ragtag
- gap2seq
- blast+
- pymlst
- augur
- ska2

### Output

| File Name                           | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        |
| ----------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| status/                             | empty .txt files produced at the end of each rule, used by snakemake to determine which rules need to be run and what order to run them in                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         |
| **For Denovo Assembly**             |
| denovo.blast.tsv                    | results of blast nucleotide search on denovo assembled sequences in ompA database to perform ompA genotyping                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
| denovo.mlst.ct.results.tsv          | genotyping of all denovo assembled sequences in chlamydiales order of PubMLST database                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             |
| denovo.mlst.generic.results.tsv     | genotyping of all denovo assembled sequences in C. trachomatis species of PubMLST database                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         |
| denovo.mlst.plasmid.results.tsv     | genotyping of all denovo assembled sequences in PubMLST's plasmid database                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         |
| ct.tree                             | phylogenetic tree of all denovo assembled sequences and reference sequences                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        |
| tree/                               | input.list: list of all denovo assembled and reference sequences along with their filepaths<br>iqtree.log: log output of tree generation<br>ska_alignment-delim.fasta.contree: phylogenetic tree topology of aligned sequences<br>ska_alignment-delim.fasta.log: log output of alignment process<br>ska_alignment-delim.fasta.splits.nex: NEXUS file format of phylogenetic tree<br>ska_alignment-delim.iqtree.log: duplicate of ska_alignment-delim.fasta.log<br>ska_alignment.fasta: all sequences aligned to reference sequences<br>ska_alignment.fasta.csv: count of each base type in every sequence (assembled and reference)<br>ska.log: log output of ska2 |
| **For Reference Assembly**          |
| ref-denovo.ompA.blast.tsv                | results of blast nucleotide search on reference assembled sequences in ompA database to perform ompA genotyping                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    |
| ref-denovo.coverage.tsv             | coverage statistics of all reference assembled sequences compared to all reference sequences                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
| ref-denovo.mlst.ct.results.tsv      | genotyping of all reference assembled sequences in chlamydiales order of PubMLST database                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |
| ref-denovo.mlst.generic.results.tsv | genotyping of all reference assembled sequences in C. trachomatis species of PubMLST database                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      |
| ref-denovo.mlst.plasmid.results.tsv | genotyping of all reference assembled sequences in PubMLST's plasmid database                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      |

**For Each Sample:**
| File Name                                                           | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
| ------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| benchmark/                                                          | information about how long each rule took                                                                                                                                                                                                                                                                                                                                                                                                                                                     |
| log/                                                                | log outputs of each step, useful for debugging                                                                                                                                                                                                                                                                                                                                                                                                                                                |
| scrub/                                                              | scrubby_temp/ kraken2 database and host genome sequences used by scrubby<br>_scrub_first_\*.fastq.gz first scrub, samples with reads classified as host, Archaea, Eukaryota, Holozoa, and Nucletmycea removed<br>_trim_scrub_\*.fastq.gz second scrub, extracted Chlamydiales reads from first scrub results<br>ct_extract.json scrub statistics after removing Archaea, Eukaryota, Holozoa, and Nucletmycea reads<br>_remove_host.json scrub statistics after removing host reads  |
| trim/                                                               | .fastq.gz trimmed samples<br>.html fastp trim statistics (number of reads and bases before and after trimming, sequence content, etc)<br>.json json version of html                                                                                                                                                                                                                                                                                                                           |
| **With denovo mode enabled, denovo/**                               |
| denovo_{sample}.final.fasta                                         | final denovo assembled sequence                                                                                                                                                                                                                                                                                                                                                                                                                                                               |
| blast/blast.ompa.tab                                                | tab separated blast nucleotide search for this particular sample (1 row of denovo.blast.tsv)                                                                                                                                                                                                                                                                                                                                                                                                  |
| gap2seq/filled.fasta                                                | filled in scaffold according to trimmed and scrubbed sequences                                                                                                                                                                                                                                                                                                                                                                                                                                |
| mlst/                                                               | {sample}.genome.chlamydiales.mlst.txt genotyping in PubMLST database of assembled sequence in chlamydiales order (1 row in denovo.mlst.generic.results.tsv)<br>{sample}.genome.ctrachomatis.mlst.txt genotyping in PubMLST database of assembled sequence in C. trachomatis species (1 row in denovo.mlst.ct.results.tsv)<br>{sample}.genome.plasmid.mlst.txt genotyping in PubMLST database of assembled sequence in plasmid database (1 row in denovo.mlst.plasmid.results.tsv)             |
| quast_report/                                                       | See [Quast Github page](https://github.com/ablab/quast?tab=readme-ov-file#output)                                                                                                                                                                                                                                                                                                                                                                                                             |
| scaffold/                                                           | See [RagTag wiki](https://github.com/malonge/RagTag/wiki/scaffold#output)                                                                                                                                                                                                                                                                                                                                                                                                                     |
| shovill/                                                            | See [Shovill Github page](https://github.com/tseemann/shovill?tab=readme-ov-file#output-files)                                                                                                                                                                                                                                                                                                                                                                                                |
| **With reference assembly mode enabled, ref-denovo/**               |
| ref-denovo_{sample}.final.fasta                                     | final reference assembled sequence                                                                                                                                                                                                                                                                                                                                                                                                                                                            |
| blast/blast.ompa.tab                                                | tab separated blast nucleotide search for this particular sample (1 row of ref-denovo.blast.tsv)                                                                                                                                                                                                                                                                                                                                                                                              |
| bowtie2/6276.bam                                                    | alignment of cleaned reads to reference sequence                                                                                                                                                                                                                                                                                                                                                                                                                                              |
| bowtie_ref24/                                                       | {sample}.coverage.ref24.tsv coverage statistics of this sample and reference sequences (part of ref-denovo.coverage.tsv)<br>{sample}.ref24.bam alignment of cleaned reads to all reference sequences<br>{sample}.ref24.bam.bai indexes for {sample}.ref24.bam<br>{sample}.tmp.cov.tsv {sample}.coverage.ref24.tsv without sample_name column                                                                                                                                                  |
| gap2seq/filled.fasta                                                | filled in scaffold according to trimmed and scrubbed sequences                                                                                                                                                                                                                                                                                                                                                                                                                                |
| mlst/                                                               | {sample}.genome.chlamydiales.mlst.txt genotyping in PubMLST database of assembled sequence in chlamydiales order (1 row in ref-denovo.mlst.generic.results.tsv)<br>{sample}.genome.ctrachomatis.mlst.txt genotyping in PubMLST database of assembled sequence in C. trachomatis species (1 row in ref-denovo.mlst.ct.results.tsv)<br>{sample}.genome.plasmid.mlst.txt genotyping in PubMLST database of assembled sequence in plasmid database (1 row in ref-denovo.mlst.plasmid.results.tsv) |
| quast_report/                                                       | See [Quast Github page](https://github.com/ablab/quast?tab=readme-ov-file#output)                                                                                                                                                                                                                                                                                                                                                                                                             |
| scaffold/                                                           | See [RagTag wiki](https://github.com/malonge/RagTag/wiki/scaffold#output)                                                                                                                                                                                                                                                                                                                                                                                                                     |
| shovill/                                                            | See [Shovill Github page](https://github.com/tseemann/shovill?tab=readme-ov-file#output-files)                                                                                                                                                                                                                                                                                                                                                                                                |
| **Additional result with both enabled**                             |
| dnadiff_report/                                                     | See [DNADiff Github page](https://github.com/garviz/MUMmer/blob/master/docs/dnadiff.README)                                                                                                                                                                                                                                                                                                                                                                                                   |

### Cite

Pipeline is created by Shola Olagoke with assistance from Ammar Aziz and Lucile Zhu.
