# ============================================================
# Makefile: RNA-Seq pipeline for UHR vs HBR chr22 subset
# ============================================================

SHELL := /usr/bin/env bash
.SHELLFLAGS := -euo pipefail -c

# -----------------------
# Config (overridable)
# -----------------------

THREADS ?= 4
DESIGN  ?= design.csv

# chr22 subset from the Biostar uhr-hbr.tar.gz package
REF := refs/chr22.genome.fa
GTF := refs/chr22.gtf

# -----------------------
# Directories
# -----------------------

BAM_DIR   := bam
BW_DIR    := bw
COUNT_DIR := counts

# -----------------------
# Samples (from design.csv)
# -----------------------

SAMPLES := $(shell awk -F, 'NR>1 && NF {print $$1}' $(DESIGN))

BAMS := $(patsubst %,$(BAM_DIR)/%.bam,$(SAMPLES))
BWS  := $(patsubst %,$(BW_DIR)/%.bw,$(SAMPLES))

# -----------------------
# Tools
# -----------------------

HISAT2        ?= hisat2
SAMTOOLS      ?= samtools
BEDTOOLS      ?= bedtools
BG2BW         ?= bedGraphToBigWig
FEATURECOUNTS ?= featureCounts

# -----------------------
# Phony targets
# -----------------------

.PHONY: all toolcheck align bigwig counts clean

# Full pipeline
all: toolcheck align bigwig counts

# -----------------------
# Tool check
# -----------------------

toolcheck:
	@echo "Checking required tools..."
	@missing=0; \
	for t in "$(HISAT2)" "$(SAMTOOLS)" "$(BEDTOOLS)" "$(BG2BW)" "$(FEATURECOUNTS)"; do \
	  if ! command -v $$t >/dev/null 2>&1; then \
	    echo "  MISSING: $$t"; missing=1; \
	  else \
	    echo "  OK: $$t"; \
	  fi; \
	done; \
	if [ $$missing -ne 0 ]; then \
	  echo "ERROR: install missing tools in the 'bioinfo' env and re-run."; \
	  exit 1; \
	fi

# -----------------------
# Directories
# -----------------------

$(BAM_DIR):
	mkdir -p $(BAM_DIR)

$(BW_DIR):
	mkdir -p $(BW_DIR)

$(COUNT_DIR):
	mkdir -p $(COUNT_DIR)

# -----------------------
# Alignment: HISAT2 -> sorted BAM + index
# -----------------------
# FASTQs are already present from uhr-hbr.tar.gz:
#   reads/HBR_1_R1.fq, reads/UHR_3_R1.fq, etc.

# pattern: sample name is e.g. HBR_1, FASTQ = reads/HBR_1_R1.fq

$(BAM_DIR)/%.bam: $(REF) reads/%_R1.fq | $(BAM_DIR)
	@echo "Aligning sample $* with HISAT2..."
	$(HISAT2) -p $(THREADS) \
	  -x $(REF) \
	  -U reads/$*_R1.fq \
	  | $(SAMTOOLS) sort -@ $(THREADS) -o $@
	$(SAMTOOLS) index $@

align: $(BAMS)

# -----------------------
# Coverage: BAM -> BigWig
# -----------------------

$(BW_DIR)/%.bw: $(BAM_DIR)/%.bam $(REF).fai | $(BW_DIR)
	@echo "Generating BigWig for $* ..."
	$(BEDTOOLS) genomecov -ibam $< -bg -split \
	  | sort -k1,1 -k2,2n > $(BW_DIR)/$*.bedgraph
	$(BG2BW) $(BW_DIR)/$*.bedgraph $(REF).fai $@
	rm -f $(BW_DIR)/$*.bedgraph

bigwig: $(BWS)

# -----------------------
# FeatureCounts: gene-level count matrix
# -----------------------

$(COUNT_DIR)/counts.txt: $(BAMS) $(GTF) | $(COUNT_DIR)
	@echo "Running featureCounts on all BAMs..."
	$(FEATURECOUNTS) -T $(THREADS) \
	  -a $(GTF) -t exon -g gene_id \
	  -o $@ $(BAMS)
	@echo "Counts written to $@"

counts: $(COUNT_DIR)/counts.txt

# -----------------------
# Cleanup
# -----------------------

clean:
	rm -rf $(BAM_DIR) $(BW_DIR) $(COUNT_DIR)
	# refs/ and reads/ are left intact on purpose