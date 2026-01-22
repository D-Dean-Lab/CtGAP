# ============================================
#  CtGAP - Chlamydia Genome Assembly Pipeline
# ============================================
.PHONY: all setup run clean help download-databases

# Platform detection - use x86 emulation on Mac ARM for bioconda compatibility
UNAME_S := $(shell uname -s)
UNAME_M := $(shell uname -m)
ifeq ($(UNAME_S),Darwin)
    # macOS - get core count and use half (leave resources for other apps)
    TOTAL_CORES := $(shell sysctl -n hw.ncpu)
    JOBS := $(shell echo $$(( $(shell sysctl -n hw.ncpu) / 2 )))
    ifeq ($(UNAME_M),arm64)
        # Mac ARM (M1/M2/M3) - use Rosetta for bioconda packages
        export CONDA_SUBDIR=osx-64
        PLATFORM_MSG := "Platform: macOS ARM64 (using x86 emulation for compatibility)"
    else
        PLATFORM_MSG := "Platform: macOS x86_64"
    endif
else
    # Linux - cap at 32 jobs to prevent memory exhaustion during assembly
    TOTAL_CORES := $(shell nproc)
    JOBS := 32
    PLATFORM_MSG := "Platform: Linux $(UNAME_M)"
endif

# Ensure at least 2 jobs
ifeq ($(shell test $(JOBS) -lt 2 && echo yes),yes)
    JOBS := 2
endif

# Config overrides from command line (via ctgap wrapper)
CONFIG_OVERRIDES :=
ifdef CTGAP_MODE
    CONFIG_OVERRIDES += mode=$(CTGAP_MODE)
endif
ifdef CTGAP_REFERENCE
    CONFIG_OVERRIDES += reference=$(CTGAP_REFERENCE)
endif

# Convert to snakemake --config format
SNAKEMAKE_CONFIG := $(if $(CONFIG_OVERRIDES),--config $(CONFIG_OVERRIDES),)

# Database URLs
GRCH38_URL := https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
KRAKEN2_URL := https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08_GB_20251015.tar.gz
BAKTA_DB_URL := https://zenodo.org/records/14025623/files/db-light.tar.gz

# Default action: full pipeline
all: setup run

# --- Database downloads ---
download-databases:
	@echo ""
	@echo "📦 [CtGAP] Checking required databases..."
	@echo ""
	@mkdir -p resources
	@# Check and download GRCh38 human reference
	@if [ -f "resources/grch38.fasta" ]; then \
		echo "   ✓ GRCh38 reference found (resources/grch38.fasta)"; \
	else \
		echo "   ⬇ Downloading GRCh38 human reference genome (~900 MB compressed)..."; \
		echo "     This may take several minutes depending on your connection."; \
		curl -L --progress-bar -o resources/grch38.fasta.gz "$(GRCH38_URL)" && \
		echo "   ⬇ Decompressing GRCh38 (~3 GB uncompressed)..." && \
		gunzip -f resources/grch38.fasta.gz && \
		echo "   ✓ GRCh38 reference downloaded successfully"; \
	fi
	@echo ""
	@# Check and download Kraken2 database
	@if [ -d "resources/standardDB" ] && [ -f "resources/standardDB/hash.k2d" ]; then \
		echo "   ✓ Kraken2 database found (resources/standardDB/)"; \
	else \
		echo "   ⬇ Downloading Kraken2 standard database (~8 GB)..."; \
		echo "     This may take 10-30 minutes depending on your connection."; \
		mkdir -p resources/standardDB && \
		curl -L --progress-bar -o resources/standardDB/k2_standard_08gb.tar.gz "$(KRAKEN2_URL)" && \
		echo "   ⬇ Extracting Kraken2 database..." && \
		tar -xzf resources/standardDB/k2_standard_08gb.tar.gz -C resources/standardDB && \
		rm -f resources/standardDB/k2_standard_08gb.tar.gz && \
		echo "   ✓ Kraken2 database downloaded successfully"; \
	fi
	@echo ""
	@echo "   ✓ Core databases ready (Bakta DB downloaded after conda setup)"
	@echo ""

# --- Environment setup ---
setup: download-databases
	@echo "🔧 [CtGAP] Setting up conda environments (first-time only)..."
	@echo $(PLATFORM_MSG)
ifdef CTGAP_MODE
	@echo "Mode: $(CTGAP_MODE)"
endif
ifdef CTGAP_REFERENCE
	@echo "Reference: $(CTGAP_REFERENCE)"
endif
	# Ensure standard folders exist
	@mkdir -p input output/status output/log output/benchmark
	# Create tiny dummy FASTQs only if no real samples exist
	@if [ -z "$$(ls input/*_R1.fastq.gz 2>/dev/null)" ]; then \
		echo "No samples found, creating DUMMY files for environment setup..."; \
		printf "@r1\nACGT\n+\n!!!!\n" | gzip -c > input/DUMMY_R1.fastq.gz; \
		printf "@r2\nTGCA\n+\n!!!!\n" | gzip -c > input/DUMMY_R2.fastq.gz; \
	fi
	# Build only the conda envs (no real compute)
	@CTGAP_ENVONLY=1 CONDA_SOLVER=classic snakemake \
		--snakefile workflow/snakefile \
		--configfile config/config.yaml \
		$(SNAKEMAKE_CONFIG) \
		--use-conda \
		--conda-create-envs-only \
		-j1 \
		--conda-frontend conda
	# Remove DUMMY files after environment setup
	@rm -f input/DUMMY_R1.fastq.gz input/DUMMY_R2.fastq.gz
	# Download Bakta database (if not present)
# Uses direct wget download (bakta_db download has checksum issues)
# URL: Bakta v6 light database from Zenodo
	@if [ -d "resources/bakta_db/db-light" ] && [ -f "resources/bakta_db/db-light/version.json" ]; then \
		echo "   ✓ Bakta database found (resources/bakta_db/db-light)"; \
	else \
		echo ""; \
		echo "   ⬇ Downloading Bakta light database (~1.3 GB)..."; \
		echo "     This may take several minutes depending on your connection."; \
		mkdir -p resources/bakta_db && \
		cd resources/bakta_db && \
		rm -f db-light.tar.xz && \
		wget --tries=3 --timeout=120 --show-progress -O db-light.tar.xz \
			"https://zenodo.org/records/14916843/files/db-light.tar.xz" && \
		echo "   ⬇ Extracting database (this may take a few minutes)..." && \
		tar -xJf db-light.tar.xz && \
		rm -f db-light.tar.xz && \
		cd ../.. ; \
		if [ -f "resources/bakta_db/db-light/version.json" ]; then \
			echo "   ✓ Bakta database downloaded successfully"; \
		else \
			echo "   ❌ Bakta database download failed!"; \
			echo "      Please check your internet connection and try again."; \
			exit 1; \
		fi; \
	fi
	# Update AMRFinderPlus database (required for Bakta annotation) - only if not present
	@if [ -L "resources/bakta_db/db-light/amrfinderplus-db/latest" ]; then \
		echo "   ✓ AMRFinderPlus database found (resources/bakta_db/db-light/amrfinderplus-db)"; \
	else \
		echo "   ⬇ Setting up AMRFinderPlus database..."; \
		BAKTA_ENV=$$(find .snakemake/conda -maxdepth 2 -type d 2>/dev/null | while read d; do \
			if [ -f "$$d/bin/amrfinder_update" ]; then echo "$$d"; break; fi; \
		done); \
		if [ -n "$$BAKTA_ENV" ]; then \
			conda run -p "$$BAKTA_ENV" amrfinder_update --force_update \
				--database resources/bakta_db/db-light/amrfinderplus-db 2>&1 | tail -3 && \
			echo "   ✓ AMRFinderPlus database setup complete"; \
		else \
			echo "   ⚠ Could not find bakta environment for AMRFinderPlus update"; \
		fi; \
	fi
	@echo ""
	@echo "✅ ═══════════════════════════════════════════════"
	@echo "   CtGAP setup complete!"
	@echo "   - Databases: ready"
	@echo "   - Conda environments: created"
	@echo "   Ready to process your samples."
	@echo "═══════════════════════════════════════════════ ✅"
	@echo ""

# --- Pipeline run ---
run:
	@echo ""
	@echo "═════════════════════════════════════════════"
	@echo "   Starting CtGAP Pipeline..."
	@echo "═════════════════════════════════════════════"
	@echo $(PLATFORM_MSG)
	@echo "Using $(JOBS) parallel jobs ($(TOTAL_CORES) cores detected)"
ifdef CTGAP_MODE
	@echo "Mode: $(CTGAP_MODE)"
endif
ifdef CTGAP_REFERENCE
	@echo "Reference: $(CTGAP_REFERENCE)"
endif
	@echo ""
	@CTGAP_ENVONLY=0 snakemake \
		--snakefile workflow/snakefile \
		--configfile config/config.yaml \
		$(SNAKEMAKE_CONFIG) \
		--use-conda \
		-j $(JOBS) \
		--rerun-incomplete \
		-q progress \
		--show-failed-logs \
		-k \
		--latency-wait 180 && \
	{ \
		echo ""; \
		echo "═══════════════════════════════════════════════"; \
		echo ""; \
		echo "  ✨ Pipeline completed successfully!"; \
		echo "  📁 Check your results in: output/"; \
		echo ""; \
		echo "═══════════════════════════════════════════════"; \
		echo ""; \
	} || \
	{ \
		echo ""; \
		echo "❌ ═══════════════════════════════════════════════"; \
		echo "   CtGAP encountered an error!"; \
		echo "   Check the logs in .snakemake/log/"; \
		echo "═══════════════════════════════════════════════ ❌"; \
		echo ""; \
		exit 1; \
	}

# --- Clean up ---
# Removes intermediate files and consolidates output to 01-04 directories
# Run this when you're satisfied with results and want to save disk space
clean:
	@echo ""
	@echo "🧹 [CtGAP] Cleaning up and finalizing outputs..."
	@echo "   This will:"
	@echo "   - Remove large intermediate files (FASTQs, BAMs, temp files)"
	@echo "   - Convert organized directories to standalone (remove symlinks)"
	@echo "   - Keep only 01-04 directories with final results"
	@echo ""
	@# Step 1: Remove snakemake cache (but preserve conda environments!)
	@rm -rf .snakemake/log .snakemake/locks .snakemake/incomplete .snakemake/shadow .snakemake/metadata .snakemake/auxiliary 2>/dev/null || true
	@# Keep .snakemake/conda/ - these take a long time to build
	@rm -f input/DUMMY_R1.fastq.gz input/DUMMY_R2.fastq.gz
	@echo "   ✓ Removed snakemake cache"
	@# Step 2: Remove large intermediate files from sample directories
	@# Trimmed and scrubbed reads
	@find output -path "*/trim/*.fastq.gz" -type f -delete 2>/dev/null || true
	@find output -path "*/trim/*.fq.gz" -type f -delete 2>/dev/null || true
	@find output -path "*/scrub/*.fastq.gz" -type f -delete 2>/dev/null || true
	@find output -path "*/scrub/*.fq.gz" -type f -delete 2>/dev/null || true
	@# Kraken2 raw output (huge - same size as input reads). Keep .report files
	@find output -name "*.kraken2.out" -type f -delete 2>/dev/null || true
	@find output -name "*.kraken2" -type f -delete 2>/dev/null || true
	@find output -name "*.kraken2.labels" -type f -delete 2>/dev/null || true
	@# BAM files (can be 100s of MB each)
	@find output -name "*.bam" -type f -delete 2>/dev/null || true
	@find output -name "*.bam.bai" -type f -delete 2>/dev/null || true
	@# Shovill temp files (denovo and ref-denovo)
	@find output -path "*/shovill/spades" -type d -exec rm -rf {} + 2>/dev/null || true
	@find output -path "*/shovill/*.fq" -type f -delete 2>/dev/null || true
	@find output -path "*/shovill/*.fq.gz" -type f -delete 2>/dev/null || true
	@# Gap2Seq temp files
	@find output -path "*/gap2seq/*.bed" -type f -delete 2>/dev/null || true
	@# Unicycler temp files (can be very large)
	@find output -path "*/unicycler/spades_assembly" -type d -exec rm -rf {} + 2>/dev/null || true
	@find output -path "*/unicycler/*.fq" -type f -delete 2>/dev/null || true
	@find output -path "*/unicycler/*.fastq" -type f -delete 2>/dev/null || true
	@# Reference-guided bamtofastq extracted reads (50-200 MB each)
	@find output -path "*/bamtofastq/*.fastq.gz" -type f -delete 2>/dev/null || true
	@find output -path "*/bamtofastq/*.fq.gz" -type f -delete 2>/dev/null || true
	@# Benchmark files
	@rm -rf output/benchmark 2>/dev/null || true
	@find output -name "benchmark" -type d -exec rm -rf {} + 2>/dev/null || true
	@# Log files (can accumulate)
	@find output -name "log" -type d -exec rm -rf {} + 2>/dev/null || true
	@echo "   ✓ Removed intermediate files (FASTQs, BAMs, temp files)"
	@# Step 3: Convert symlinks to real files in 01-04 directories
	@echo "   ✓ Converting symlinks to standalone files..."
	@for dir in output/01_per_sample output/02_reports output/03_results output/04_phylogeny; do \
		if [ -d "$$dir" ]; then \
			find "$$dir" -type l | while read link; do \
				target=$$(readlink -f "$$link" 2>/dev/null); \
				if [ -e "$$target" ]; then \
					rm "$$link"; \
					if [ -d "$$target" ]; then \
						cp -r "$$target" "$$link"; \
					else \
						cp "$$target" "$$link"; \
					fi; \
				fi; \
			done; \
		fi; \
	done
	@echo "   ✓ Converted symlinks to real files"
	@# Step 4: Remove original root-level files/directories (now duplicated in 01-04)
	@# Keep: 00_README.txt, 01-04 directories
	@# Remove: sample dirs, root TSVs, reports, tree, status, logs, ct_typing, etc.
	@for item in output/*/; do \
		dirname=$$(basename "$$item"); \
		case "$$dirname" in \
			00_README.txt|01_per_sample|02_reports|03_results|04_phylogeny) ;; \
			*) rm -rf "$$item" 2>/dev/null || true ;; \
		esac; \
	done
	@rm -f output/*.tsv output/*.tree 2>/dev/null || true
	@rm -rf output/status output/log output/logs output/tree output/ct_typing output/multiqc 2>/dev/null || true
	@echo "   ✓ Removed original files (data now in 01-04 directories)"
	@echo ""
	@echo "✅ Cleanup complete!"
	@echo ""
	@echo "   Final output structure:"
	@ls -la output/ 2>/dev/null | grep -E "^d|README" | head -10 || true
	@echo ""
	@du -sh output 2>/dev/null || true
	@echo ""
	@echo "   ⚠️  Note: To re-run the pipeline, you'll need to start from scratch."

# --- Help message ---
help:
	@echo ""
	@echo "════════════════════════════════════════════════════"
	@echo "   CtGAP - Chlamydia trachomatis Genome Assembly"
	@echo "════════════════════════════════════════════════════"
	@echo ""
	@echo "  Commands:"
	@echo "    ctgap setup              - Download databases + create conda environments"
	@echo "    ctgap run                - Run the assembly pipeline"
	@echo "    ctgap clean              - Finalize outputs and remove intermediate files"
	@echo "    ctgap help               - Show this help message"
	@echo ""
	@echo "  Setup Details:"
	@echo "    'ctgap setup' will automatically download (if not present):"
	@echo "    - GRCh38 human reference genome (~3 GB)"
	@echo "    - Kraken2 standard database (~8 GB)"
	@echo "    - Bakta light database (~1.5 GB) for genome annotation"
	@echo "    - All required conda environments"
	@echo ""
	@echo "  Assembly Mode Options:"
	@echo "    ctgap run                            # Auto mode (default)"
	@echo "    ctgap run --denovo                   # De novo assembly only"
	@echo "    ctgap run --ref-guided E_bour        # Reference-guided with E_bour"
	@echo "    ctgap run --ref-guided plurality     # Reference-guided with plurality"
	@echo ""
	@echo "  Cleanup (run after reviewing results):"
	@echo "    ctgap clean - Removes intermediate files and consolidates output"
	@echo "                  to 01-04 directories only. Saves ~80-90% disk space."
	@echo "                  WARNING: Re-running pipeline requires starting fresh."
	@echo ""
	@echo "  Output Structure (after pipeline completes):"
	@echo "    output/01_per_sample/  - Per-sample data"
	@echo "    output/02_reports/     - All reports (PDFs, summaries)"
	@echo "    output/03_results/     - Collated typing/coverage results"
	@echo "    output/04_phylogeny/   - Phylogenetic tree"
	@echo ""
	@echo "════════════════════════════════════════════════════"
	@echo ""

# ============================================
# Developer/Maintainer targets (commented out)
# ============================================
# Uncomment these if you need to regenerate pipeline workflow/visualizations for documentation

# dag:
# 	@echo "🔍 [CtGAP] Generating workflow DAG..."
# 	snakemake \
# 		--snakefile workflow/snakefile \
# 		--configfile config/config.yaml \
# 		--dag | dot -Tpdf > dag.pdf
# 	@echo "✅ Full DAG saved to dag.pdf"
