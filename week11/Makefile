# Use bash and fail fast
SHELL := /usr/bin/env bash
.SHELLFLAGS := -euo pipefail -c

# -------------------------------------
# Variables (overridable)
# -------------------------------------
REF_ASM ?= GCF_000848505.1
THREADS ?= 4
JOBS ?= 4
SUBSET ?=
DESIGN ?= design.csv

# Reproducibility guard:
# STRICT=1 -> require EVERY sample in design.csv to produce a VCF, else fail
# STRICT=0 -> merge only VCFs that exist (debug mode)
STRICT ?= 1

# Column names in design.csv
DESIGN_RUN ?= Run
DESIGN_SAMPLE ?= Sample
DESIGN_LAYOUT ?= Layout

# Per-sample variables
SRR ?= SRR1734993
SAMPLE ?= $(SRR)
LAYOUT ?= PE

# Variant calling
PLOIDY ?= 1
MIN_MAPQ ?= 20
MIN_BASEQ ?= 20

# Genome download
DOWNLOAD_RETRIES ?= 3
SLEEP_BETWEEN ?= 5
CURL ?= curl -L
ESEARCH ?= esearch
XTRACT ?= xtract
GUNZIP ?= gunzip

# -------------------------------------
# Directories
# -------------------------------------
GENOME_DIR = genome
READS_DIR = reads
RESULTS_DIR = results/$(SAMPLE)
MERGED_DIR = results/merged
JOBS_DIR   = jobs
SRA_CACHE  = $(HOME)/ncbi/public/sra

# -------------------------------------
# Files
# -------------------------------------
GENOME = $(GENOME_DIR)/ebola_genome.fna
GFF = $(GENOME_DIR)/ebola.gff

BAM = $(RESULTS_DIR)/$(SAMPLE)_sorted.bam
BAI = $(BAM).bai
STATS = $(RESULTS_DIR)/$(SAMPLE)_alignment_stats.txt

BW = $(RESULTS_DIR)/$(SAMPLE)_coverage.bw
BG = $(RESULTS_DIR)/$(SAMPLE)_coverage.bedgraph
CHROMSIZES = $(GENOME).chrom.sizes

READ1 = $(READS_DIR)/$(SRR)_1.fastq
READ2 = $(READS_DIR)/$(SRR)_2.fastq
READ_SE = $(READS_DIR)/$(SRR).fastq

VCF_GZ = $(RESULTS_DIR)/$(SAMPLE).vcf.gz
VCF_TBI = $(VCF_GZ).tbi
VCF_STATS = $(RESULTS_DIR)/$(SAMPLE)_vcf_stats.txt

MERGED_VCF = $(MERGED_DIR)/merged.vcf.gz
MERGED_TBI = $(MERGED_VCF).tbi
VCF_LIST = $(MERGED_DIR)/vcfs.list

ifeq ($(LAYOUT),PE)
  FASTQ_TARGET = $(READ1) $(READ2)
else
  FASTQ_TARGET = $(READ_SE)
endif

# -------------------------------------
# Tools
# -------------------------------------
BWA = bwa
SAMTOOLS = samtools
BEDTOOLS = bedtools
BG2BW = bedGraphToBigWig
FASTQC = fastqc
PREFETCH = prefetch
FASTQ_DUMP = fastq-dump
FASTERQ_DUMP = fasterq-dump
DATASETS = datasets
UNZIP = unzip
BCFTOOLS = bcftools
TABIX = tabix
PARALLEL ?= parallel --eta --bar --lb --keep-order

# --- snpEff-only annotation ---
SNPEFF ?= snpEff
SNPEFF_DATA ?= snpeff_data
SNPEFF_GENOME_ID ?= ebola
SNPEFF_CONFIG ?= $(SNPEFF_DATA)/snpEff.config

VCF_SNPEFF = $(RESULTS_DIR)/$(SAMPLE).snpeff.vcf.gz
VCF_SNPEFF_TBI = $(VCF_SNPEFF).tbi
EFFECTS_SUMMARY = $(RESULTS_DIR)/$(SAMPLE)_effects_summary.txt
MERGED_SNPEFF = $(MERGED_DIR)/merged.snpeff.vcf.gz

.PHONY: help run setup-genome index fastq qc align stats bigwig vcf vcf_stats all batch one clean toolcheck \
        multisample vcf_list vcf_list_strict _ensure_dirs design_check jobsdir annotate effects snpeff_db effects_snpeff effects_summary annotate_merged manifest

# Choose list builder based on STRICT
ifeq ($(STRICT),1)
  VCF_LIST_TARGET = vcf_list_strict
else
  VCF_LIST_TARGET = vcf_list
endif

# -------------------------------------
# Help
# -------------------------------------
run: batch
	@echo "Done. Outputs under results/<Sample>/"

help:
	@echo "Targets:"
	@echo "  batch           - Validate design, run all samples (Parallel), merge VCFs"
	@echo "  one             - Run single sample (SRR=..., SAMPLE=..., LAYOUT=PE|SE)"
	@echo "  annotate        - snpEff annotate per-sample VCF"
	@echo "  effects         - snpEff annotate + summarize >=3 effect types"
	@echo "  annotate_merged - snpEff annotate merged VCF"
	@echo "  STRICT=1        - require all per-sample VCFs (default)"
	@echo "  STRICT=0        - merge only existing VCFs (debug mode)"
	@echo "  clean           - Remove reads/, results/, genome/, jobs/, snpeff_data/"

# -------------------------------------
# Tool check
# -------------------------------------
toolcheck:
	echo "Checking required tools in PATH..."
	missing=0; \
	for t in "$(firstword $(PARALLEL))" $(BWA) $(SAMTOOLS) $(BEDTOOLS) $(BG2BW) $(FASTQC) $(PREFETCH) $(FASTQ_DUMP) $(DATASETS) $(UNZIP) $(BCFTOOLS) $(TABIX) $(SNPEFF); do \
	  if ! command -v $$t >/dev/null 2>&1; then echo "  MISSING: $$t"; missing=1; else echo "  OK: $$t"; fi; \
	done; \
	if [ $$missing -ne 0 ]; then echo "ERROR: Install missing tools."; exit 127; fi

# -------------------------------------
# Genome download
# -------------------------------------
$(GENOME):
	@echo "Downloading reference genome for $(REF_ASM)..."
	mkdir -p $(GENOME_DIR)
	ok=0; \
	for i in $$(seq 1 $(DOWNLOAD_RETRIES)); do \
	  echo "Attempt $$i via NCBI datasets..."; \
	  if $(DATASETS) download genome accession $(REF_ASM) --include genome,gff3 --filename $(GENOME_DIR)/$(REF_ASM).zip; then ok=1; break; fi; \
	  echo "datasets failed (attempt $$i). Retrying..."; \
	  sleep $(SLEEP_BETWEEN); \
	done; \
	if [ $$ok -eq 1 ]; then \
	  $(UNZIP) -o $(GENOME_DIR)/$(REF_ASM).zip -d $(GENOME_DIR)/tmp_unzip; \
	  mv $(GENOME_DIR)/tmp_unzip/ncbi_dataset/data/*/*.fna $(GENOME); \
	  if compgen -G "$(GENOME_DIR)/tmp_unzip/ncbi_dataset/data/*/*.gff" > /dev/null; then mv $(GENOME_DIR)/tmp_unzip/ncbi_dataset/data/*/*.gff $(GFF); fi; \
	  rm -rf $(GENOME_DIR)/tmp_unzip $(GENOME_DIR)/$(REF_ASM).zip; \
	else \
	  echo "Failed to download genome after $(DOWNLOAD_RETRIES) attempts."; exit 1; \
	fi
	@echo "Genome written to: $(GENOME)"

setup-genome: $(GENOME)
	@echo "Indexing genome..."
	$(BWA) index $(GENOME)
	$(SAMTOOLS) faidx $(GENOME)
	@echo "Genome indexing complete."

# -------------------------------------
# Download reads
# -------------------------------------
fastq: $(FASTQ_TARGET)

$(FASTQ_TARGET): _ensure_dirs
	@[ -n "$(SRR)" ] || { echo "ERROR: SRR is empty"; exit 64; }
	@[ -n "$(SAMPLE)" ] || { echo "ERROR: SAMPLE is empty"; exit 64; }
	@[ -n "$(LAYOUT)" ] || { echo "ERROR: LAYOUT is empty"; exit 64; }
	@echo "Downloading reads for $(SRR) (LAYOUT=$(LAYOUT))..."
	mkdir -p $(READS_DIR)
	$(PREFETCH) -O $(SRA_CACHE) $(SRR) || true
	if [ -s "$(SRA_CACHE)/$(SRR).sra" ]; then SRC="$(SRA_CACHE)/$(SRR).sra"; else SRC="$(SRR)"; fi; \
	if [ "$(LAYOUT)" = "PE" ]; then \
	  $(FASTQ_DUMP) --split-files -O $(READS_DIR) "$$SRC" || \
	  $(FASTERQ_DUMP) --split-files -S -e $(THREADS) -O $(READS_DIR) "$$SRC"; \
	else \
	  $(FASTQ_DUMP) -O $(READS_DIR) "$$SRC" || \
	  $(FASTERQ_DUMP) -S -e $(THREADS) -O $(READS_DIR) "$$SRC"; \
	fi

# -------------------------------------
# QC
# -------------------------------------
qc: fastq
	@echo "Running FastQC for $(SRR)..."
	inputs=""; \
	if [ -f "$(READ1)" ]; then inputs="$$inputs $(READ1)"; fi; \
	if [ -f "$(READ2)" ]; then inputs="$$inputs $(READ2)"; fi; \
	if [ -f "$(READ_SE)" ]; then inputs="$$inputs $(READ_SE)"; fi; \
	if [ -n "$$inputs" ]; then $(FASTQC) -t $(THREADS) -o $(READS_DIR) $$inputs; else echo "No FASTQs found"; exit 1; fi

# -------------------------------------
# Alignment
# -------------------------------------
align: index fastq
	@echo "Aligning $(SRR)..."
	mkdir -p $(RESULTS_DIR)
	if [ "$(LAYOUT)" = "PE" ]; then \
	  $(BWA) mem -t $(THREADS) $(GENOME) $(READ1) $(READ2) | $(SAMTOOLS) sort -o $(BAM); \
	else \
	  $(BWA) mem -t $(THREADS) $(GENOME) $(READ_SE) | $(SAMTOOLS) sort -o $(BAM); \
	fi
	$(SAMTOOLS) index $(BAM)

# -------------------------------------
# Stats
# -------------------------------------
stats: align
	$(SAMTOOLS) flagstat $(BAM) > $(STATS)
	$(SAMTOOLS) depth $(BAM) | awk '{s+=$$3; n++} END{if(n>0) printf "Average coverage = %.3f\n", s/n; else print "Average coverage = 0"}' >> $(STATS)

# -------------------------------------
# BigWig (skips empty BAM)
# -------------------------------------
bigwig: align
	if [ "$$($(SAMTOOLS) view -c $(BAM))" -eq 0 ]; then \
	  echo "Skipping BigWig for $(SAMPLE) (no reads)"; \
	else \
	  $(SAMTOOLS) faidx $(GENOME); \
	  cut -f1,2 $(GENOME).fai > $(CHROMSIZES); \
	  $(BEDTOOLS) genomecov -ibam $(BAM) -bg | sort -k1,1 -k2,2n > $(BG); \
	  $(BG2BW) $(BG) $(CHROMSIZES) $(BW); \
	fi

# -------------------------------------
# Variant calling
# -------------------------------------
vcf: align
	@echo "Calling variants for $(SAMPLE)..."
	mkdir -p $(RESULTS_DIR)
	$(BCFTOOLS) mpileup -Ou -f $(GENOME) -q $(MIN_MAPQ) -Q $(MIN_BASEQ) -a DP,AD $(BAM) | \
	$(BCFTOOLS) call -mv --ploidy $(PLOIDY) -Ou | \
	$(BCFTOOLS) norm -f $(GENOME) -Oz -o $(VCF_GZ)
	$(TABIX) -p vcf $(VCF_GZ)
	$(BCFTOOLS) stats $(VCF_GZ) > $(VCF_STATS)

# -------------------------------------
# All per-sample
# -------------------------------------
all: qc stats bigwig vcf

# -------------------------------------
# Batch (validated + logged + manifest)
# -------------------------------------
batch: toolcheck design_check manifest
	if [ ! -f "$(GENOME)" ]; then \
	  echo "Genome not found; downloading..."; \
	  $(MAKE) setup-genome REF_ASM=$(REF_ASM) THREADS=$(THREADS); \
	fi
	@echo "Running batch for all samples in $(DESIGN)..."
	$(MAKE) jobsdir
	{ \
	  awk 'NR==1{print; next} NF' $(DESIGN) | perl -pe 's/\r//;/^\s*$$/d'; \
	} | $(PARALLEL) --colsep ',' --header : -j $(JOBS) --joblog $(JOBS_DIR)/joblog.tsv --results $(JOBS_DIR)/parallel_out \
	 '$(MAKE) all SRR={$(DESIGN_RUN)} SAMPLE={$(DESIGN_SAMPLE)} LAYOUT={$(DESIGN_LAYOUT)} THREADS=$(THREADS) SUBSET=$(SUBSET) REF_ASM=$(REF_ASM)'
	@echo "Building VCF list (STRICT=$(STRICT))..."
	$(MAKE) $(VCF_LIST_TARGET)
	$(MAKE) multisample

# Validate design.csv header contains required columns
design_check:
	@if [ ! -s "$(DESIGN)" ]; then echo "ERROR: $(DESIGN) not found or empty"; exit 64; fi
	@hdr="$$(head -n1 $(DESIGN))"; \
	case "$$hdr" in \
		*$(DESIGN_RUN)*$(DESIGN_SAMPLE)*$(DESIGN_LAYOUT)*|*$(DESIGN_RUN)*$(DESIGN_LAYOUT)*$(DESIGN_SAMPLE)*) ;; \
		*) echo "ERROR: $(DESIGN) must contain columns: $(DESIGN_RUN),$(DESIGN_SAMPLE),$(DESIGN_LAYOUT)"; echo "Found header: $$hdr"; exit 65 ;; \
	esac

jobsdir:
	mkdir -p $(JOBS_DIR) $(MERGED_DIR)

# -------------------------------------
# Build VCF list (NON-STRICT: include only existing files)
# -------------------------------------
vcf_list:
	mkdir -p $(MERGED_DIR)
	awk -F',' 'NR>1 && NF {gsub(/\r/,""); printf "results/%s/%s.vcf.gz\n",$$2,$$2}' $(DESIGN) \
	  | while read p; do [ -s "$$p" ] && echo "$$p"; done > $(VCF_LIST)
	@echo "VCFs included (existing only):"
	@cat $(VCF_LIST) || true

# -------------------------------------
# Build VCF list (STRICT: require all; fail + print exact rerun commands)
# -------------------------------------
vcf_list_strict:
	mkdir -p $(MERGED_DIR)
	awk -F',' 'NR>1 && NF {gsub(/\r/,""); printf "results/%s/%s.vcf.gz\t%s\t%s\n",$$2,$$2,$$1,$$3}' $(DESIGN) > $(MERGED_DIR)/vcfs_expected.tsv
	@missing=0; \
	> $(VCF_LIST); \
	while IFS=$$'\t' read -r vcf srr layout; do \
	  if [ -s "$$vcf" ]; then echo "$$vcf" >> $(VCF_LIST); else \
	    echo "MISSING: $$vcf  (rerun: make all SRR=$$srr SAMPLE=$$(basename $${vcf%.vcf.gz}) LAYOUT=$$layout THREADS=$(THREADS) REF_ASM=$(REF_ASM))"; \
	    missing=1; \
	  fi; \
	done < $(MERGED_DIR)/vcfs_expected.tsv; \
	if [ $$missing -ne 0 ]; then \
	  echo "ERROR: One or more per-sample VCFs missing. See MISSING lines above."; \
	  exit 66; \
	fi
	@echo "All per-sample VCFs present."

# -------------------------------------
# Merge VCFs
# -------------------------------------
multisample:
	$(BCFTOOLS) merge -l $(VCF_LIST) -Oz -o $(MERGED_VCF)
	$(TABIX) -p vcf $(MERGED_VCF)

# -------------------------------------
# snpEff-only EFFECT evaluation
# -------------------------------------
effects: annotate effects_summary
	@echo "Effects evaluation complete for $(SAMPLE). See: $(EFFECTS_SUMMARY)"

# Build a local snpEff DB from your FASTA+GFF (once)
snpeff_db: $(GENOME) $(GFF)
	@echo "Setting up snpEff DB at $(SNPEFF_DATA) for genome $(SNPEFF_GENOME_ID)..."
	mkdir -p $(SNPEFF_DATA)/data/$(SNPEFF_GENOME_ID)
	cp -f $(GENOME) $(SNPEFF_DATA)/data/$(SNPEFF_GENOME_ID)/sequences.fa
	cp -f $(GFF) $(SNPEFF_DATA)/data/$(SNPEFF_GENOME_ID)/genes.gff
	echo "# Auto-generated snpEff config" > $(SNPEFF_CONFIG)
	echo "data.dir = $(abspath $(SNPEFF_DATA))/data" >> $(SNPEFF_CONFIG)
	echo "$(SNPEFF_GENOME_ID).genome : $(SNPEFF_GENOME_ID)" >> $(SNPEFF_CONFIG)
	@echo "Building snpEff DB..."
	$(SNPEFF) build -gff3 -v -c $(SNPEFF_CONFIG) $(SNPEFF_GENOME_ID)

annotate: vcf snpeff_db
	@echo "Annotating with snpEff..."
	$(SNPEFF) -c $(SNPEFF_CONFIG) -v $(SNPEFF_GENOME_ID) $(VCF_GZ) | bgzip -c > $(VCF_SNPEFF)
	$(TABIX) -p vcf $(VCF_SNPEFF)

effects_snpeff: annotate ; @true

# Pull >=3 distinct effect types from ANN
effects_summary: $(VCF_SNPEFF)
	@echo "Summarizing effects from snpEff ANN..." > $(EFFECTS_SUMMARY)
	@echo "Input VCF: $(VCF_SNPEFF)" >> $(EFFECTS_SUMMARY)
	@echo "Field used: ANN" >> $(EFFECTS_SUMMARY)
	@echo "" >> $(EFFECTS_SUMMARY)
	$(BCFTOOLS) query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/ANN\n' $(VCF_SNPEFF) | \
	  awk -F'\t' 'function eff1(ann){split(ann,a,"|"); return a[2]} \
	              {split($$5,arr,"[,;]"); e=eff1(arr[1]); if(!(e in seen)){seen[e]=$$0}} \
	              END{printed=0; want["missense_variant"]=1; want["synonymous_variant"]=1; want["intergenic_region"]=1; \
	                  for(k in want) if(k in seen){print k"\t"seen[k]; printed++} \
	                  for(k in seen) if(printed<3 && !want[k]){print k"\t"seen[k]; printed++}}' \
	  >> $(EFFECTS_SUMMARY)
	@echo "" >> $(EFFECTS_SUMMARY)
	@echo "Notes:" >> $(EFFECTS_SUMMARY)
	@echo " - First instance per effect type shown (ANN parsed)." >> $(EFFECTS_SUMMARY)
	@echo " - For IGV, load $(GENOME), $(BAM), and $(VCF_SNPEFF)." >> $(EFFECTS_SUMMARY)

# snpEff on merged VCF (optional)
annotate_merged: multisample snpeff_db
	@echo "Annotating merged VCF with snpEff..."
	$(SNPEFF) -c $(SNPEFF_CONFIG) -v $(SNPEFF_GENOME_ID) $(MERGED_VCF) | bgzip -c > $(MERGED_SNPEFF)
	$(TABIX) -p vcf $(MERGED_SNPEFF)

# -------------------------------------
# Reproducibility manifest
# -------------------------------------
manifest:
	mkdir -p $(JOBS_DIR)
	{ \
	  echo "DATE: $$(date -u)"; \
	  echo "REF_ASM=$(REF_ASM)"; \
	  echo "THREADS=$(THREADS)  JOBS=$(JOBS)  STRICT=$(STRICT)"; \
	  echo "GENOME=$$(readlink -f $(GENOME) 2>/dev/null || echo $(GENOME))"; \
	  echo "GFF=$$(readlink -f $(GFF) 2>/dev/null || echo $(GFF))"; \
	  echo "MD5(FASTA)=$$(test -f $(GENOME) && md5sum $(GENOME) | awk '{print $$1}' || echo NA)"; \
	  echo "MD5(GFF)=$$(test -f $(GFF) && md5sum $(GFF) | awk '{print $$1}' || echo NA)"; \
	  echo "--- tool versions ---"; \
	  echo "bwa: $$(bwa 2>&1 | head -n1)"; \
	  echo "samtools: $$(samtools --version | head -n1)"; \
	  echo "bcftools: $$(bcftools --version | head -n1)"; \
	  echo "bedtools: $$(bedtools --version 2>/dev/null || echo NA)"; \
	  echo "fastqc: $$(fastqc --version 2>&1 | head -n1)"; \
	  echo "snpEff: $$(snpEff -version 2>/dev/null || echo NA)"; \
	} > $(JOBS_DIR)/manifest.txt

# -------------------------------------
# Single sample helper
# -------------------------------------
one:
	$(MAKE) all SRR=$(SRR) SAMPLE=$(SAMPLE) LAYOUT=$(LAYOUT) THREADS=$(THREADS) SUBSET=$(SUBSET) REF_ASM=$(REF_ASM)

# -------------------------------------
# Utility
# -------------------------------------
_ensure_dirs:
	mkdir -p $(READS_DIR) $(GENOME_DIR) $(MERGED_DIR)

clean:
	rm -rf $(READS_DIR) results $(GENOME_DIR) $(JOBS_DIR) $(SNPEFF_DATA)
	rm -f $(GENOME).fai $(CHROMSIZES)