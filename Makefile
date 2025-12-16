# ============================================
#  CtGAP - Chlamydia Genome Assembly Pipeline
# ============================================
.PHONY: all setup run clean help

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

# Default action: full pipeline
all: setup run

# --- Environment setup ---
setup:
	@echo "🔧 [CtGAP] Setting up environments (first-time only)..."
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
	@echo ""
	@echo "✅ ═══════════════════════════════════════════════"
	@echo "   CtGAP environments successfully created!"
	@echo "   Ready to process your samples."
	@echo "═══════════════════════════════════════════════ ✅"
	@echo ""

# --- Pipeline run ---
run:
	@echo ""
	@echo "═════════════════════════════════════════════"
	@echo "   Starting CtGAP Pipeline..."
	@echo "═════════════════════════════════════════════"
	@echo ""
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
		-j 40 \
		--rerun-incomplete \
		-p \
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
clean:
	@echo "🧹 [CtGAP] Cleaning intermediate files..."
	@rm -rf .snakemake */*.log */*/*log
	@rm -f input/DUMMY_R1.fastq.gz input/DUMMY_R2.fastq.gz
	@echo "✅ Clean complete!"

# --- Help message ---
help:
	@echo ""
	@echo "════════════════════════════════════════════════════"
	@echo "   CtGAP - Chlamydia trachomatis Genome Assembly"
	@echo "════════════════════════════════════════════════════"
	@echo ""
	@echo "  Commands:"
	@echo "    ctgap setup   - Create conda environments (first time)"
	@echo "    ctgap run     - Run the assembly pipeline"
	@echo "    ctgap clean   - Remove temporary files"
	@echo "    ctgap help    - Show this help message"
	@echo ""
	@echo "  Assembly Mode Options:"
	@echo "    ctgap run                            # Auto mode (default; runs both denovo and reference-guided, selects best for downstream)"
	@echo "    ctgap run --denovo                   # De novo assembly only"
	@echo "    ctgap run --ref-guided E_bour        # Reference guided with E_bour reference"
	@echo "    ctgap run --ref-guided plurality     # Reference-guided with plurality consensus (21 Ct strains)"
	@echo ""
	@echo "  Note: Replace 'E_bour' with any reference genome name from resources/references/ct"
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
