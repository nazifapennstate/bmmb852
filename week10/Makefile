# -------------------------------------
# Variables (overridable)
# -------------------------------------
REF_ASM      ?= GCF_000848505.1
THREADS      ?= 4
JOBS         ?= 4
SUBSET       ?=
DESIGN       ?= design.csv
DRYRUN       ?= 0

# Column names in design.csv
DESIGN_RUN    ?= Run
DESIGN_SAMPLE ?= Sample
DESIGN_LAYOUT ?= Layout

# Per-sample variables
SRR          ?= SRR1734993
SAMPLE       ?= $(SRR)
LAYOUT       ?= PE

# Variant calling
PLOIDY       ?= 1
MIN_MAPQ     ?= 20
MIN_BASEQ    ?= 20

# Genome download
DOWNLOAD_RETRIES ?= 3
SLEEP_BETWEEN    ?= 5
CURL            ?= curl -L
ESEARCH         ?= esearch
XTRACT          ?= xtract
GUNZIP          ?= gunzip

# -------------------------------------
# Directories
# -------------------------------------
GENOME_DIR   = genome
READS_DIR    = reads
RESULTS_DIR  = results/$(SAMPLE)
MERGED_DIR   = results/merged
SRA_CACHE    = $(HOME)/ncbi/public/sra

# -------------------------------------
# Files
# -------------------------------------
GENOME       = $(GENOME_DIR)/ebola_genome.fna
GFF          = $(GENOME_DIR)/ebola.gff
BAM          = $(RESULTS_DIR)/$(SAMPLE)_sorted.bam
BAI          = $(BAM).bai
STATS        = $(RESULTS_DIR)/$(SAMPLE)_alignment_stats.txt
BW           = $(RESULTS_DIR)/$(SAMPLE)_coverage.bw
BG           = $(RESULTS_DIR)/$(SAMPLE)_coverage.bedgraph
CHROMSIZES   = $(GENOME).chrom.sizes
READ1        = $(READS_DIR)/$(SRR)_1.fastq
READ2        = $(READS_DIR)/$(SRR)_2.fastq
READ_SE      = $(READS_DIR)/$(SRR).fastq
VCF_GZ       = $(RESULTS_DIR)/$(SAMPLE).vcf.gz
VCF_TBI      = $(VCF_GZ).tbi
VCF_STATS    = $(RESULTS_DIR)/$(SAMPLE)_vcf_stats.txt
MERGED_VCF   = $(MERGED_DIR)/merged.vcf.gz
MERGED_TBI   = $(MERGED_VCF).tbi
VCF_LIST     = $(MERGED_DIR)/vcfs.list

ifeq ($(LAYOUT),PE)
  FASTQ_TARGET = $(READ1) $(READ2)
else
  FASTQ_TARGET = $(READ_SE)
endif

# -------------------------------------
# Tools
# -------------------------------------
BWA          = bwa
SAMTOOLS     = samtools
BEDTOOLS     = bedtools
BG2BW        = bedGraphToBigWig
FASTQC       = fastqc
PREFETCH     = prefetch
FASTQ_DUMP   = fastq-dump
FASTERQ_DUMP = fasterq-dump
DATASETS     = datasets
UNZIP        = unzip
BCFTOOLS     = bcftools
TABIX        = tabix
PARALLEL     ?= parallel --eta --bar --lb --keep-order

.PHONY: help run setup-genome index fastq qc align stats bigwig vcf vcf_stats all batch one clean toolcheck multisample vcf_list _ensure_dirs

# -------------------------------------
# Help
# -------------------------------------
run: batch
	@echo "Done. Outputs under results/<Sample>/"

help:
	@echo "Targets:"
	@echo "  run            - Run batch over design.csv and merge results"
	@echo "  batch          - Run all rows in design.csv using GNU Parallel"
	@echo "  one            - Run single sample (SRR=..., SAMPLE=..., LAYOUT=PE|SE)"
	@echo "  multisample    - Merge all per-sample VCFs into one"
	@echo "  clean          - Remove reads/, results/, genome/"

# -------------------------------------
# Tool check
# -------------------------------------
toolcheck:
	@echo "Checking required tools in PATH..."
	@missing=0; \
	for t in "$(firstword $(PARALLEL))" $(BWA) $(SAMTOOLS) $(BEDTOOLS) $(BG2BW) $(FASTQC) $(PREFETCH) $(FASTQ_DUMP) $(DATASETS) $(UNZIP) $(BCFTOOLS) $(TABIX); do \
	  if ! command -v $$t >/dev/null 2>&1; then echo "  MISSING: $$t"; missing=1; else echo "  OK: $$t"; fi; \
	done; \
	if [ $$missing -ne 0 ]; then echo "ERROR: Install missing tools."; exit 127; fi

# -------------------------------------
# Genome download
# -------------------------------------
$(GENOME):
	@echo "Downloading reference genome for $(REF_ASM)..."
	@mkdir -p $(GENOME_DIR)
	@ok=0; \
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
	@$(BWA) index $(GENOME)
	@$(SAMTOOLS) faidx $(GENOME)
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
	@mkdir -p $(READS_DIR)
	@$(PREFETCH) -O $(SRA_CACHE) $(SRR) || true
	@if [ -s "$(SRA_CACHE)/$(SRR).sra" ]; then SRC="$(SRA_CACHE)/$(SRR).sra"; else SRC="$(SRR)"; fi; \
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
	@inputs=""; \
	if [ -f "$(READ1)" ]; then inputs="$$inputs $(READ1)"; fi; \
	if [ -f "$(READ2)" ]; then inputs="$$inputs $(READ2)"; fi; \
	if [ -f "$(READ_SE)" ]; then inputs="$$inputs $(READ_SE)"; fi; \
	if [ -n "$$inputs" ]; then $(FASTQC) -t $(THREADS) -o $(READS_DIR) $$inputs; else echo "No FASTQs found"; exit 1; fi

# -------------------------------------
# Alignment
# -------------------------------------
align: index fastq
	@echo "Aligning $(SRR)..."
	@mkdir -p $(RESULTS_DIR)
	@if [ "$(LAYOUT)" = "PE" ]; then \
	  $(BWA) mem -t $(THREADS) $(GENOME) $(READ1) $(READ2) | $(SAMTOOLS) sort -o $(BAM); \
	else \
	  $(BWA) mem -t $(THREADS) $(GENOME) $(READ_SE) | $(SAMTOOLS) sort -o $(BAM); \
	fi
	@$(SAMTOOLS) index $(BAM)

# -------------------------------------
# Stats
# -------------------------------------
stats: align
	@$(SAMTOOLS) flagstat $(BAM) > $(STATS)
	@$(SAMTOOLS) depth $(BAM) | awk '{s+=$$3; n++} END{if(n>0) printf "Average coverage = %.3f\n", s/n; else print "Average coverage = 0"}' >> $(STATS)

# -------------------------------------
# BigWig (skips empty BAM)
# -------------------------------------
bigwig: align
	@if [ "$$($(SAMTOOLS) view -c $(BAM))" -eq 0 ]; then \
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
	@mkdir -p $(RESULTS_DIR)
	@$(BCFTOOLS) mpileup -Ou -f $(GENOME) -q $(MIN_MAPQ) -Q $(MIN_BASEQ) -a DP,AD $(BAM) | \
	$(BCFTOOLS) call -mv --ploidy $(PLOIDY) -Ou | \
	$(BCFTOOLS) norm -f $(GENOME) -Oz -o $(VCF_GZ)
	@$(TABIX) -p vcf $(VCF_GZ)

vcf_stats: vcf
	@$(BCFTOOLS) stats $(VCF_GZ) > $(VCF_STATS)

# -------------------------------------
# All per-sample
# -------------------------------------
all: qc stats bigwig vcf vcf_stats

# -------------------------------------
# Batch (no temp file; stream-clean design.csv)
# -------------------------------------
batch: toolcheck
	@if [ ! -f "$(GENOME)" ]; then \
	  echo "Genome not found; downloading..."; \
	  $(MAKE) setup-genome REF_ASM=$(REF_ASM) THREADS=$(THREADS); \
	fi
	@echo "Running batch for all samples in $(DESIGN)..."
	@{ \
	  awk 'NR==1{print; next} NF' $(DESIGN) | perl -pe 's/\r//;/^\s*$$/d'; \
	} | $(PARALLEL) --colsep ',' --header : --halt now,fail=0 -j $(JOBS) \
	  '$(MAKE) all SRR={$(DESIGN_RUN)} SAMPLE={$(DESIGN_SAMPLE)} LAYOUT={$(DESIGN_LAYOUT)} THREADS=$(THREADS) SUBSET=$(SUBSET) REF_ASM=$(REF_ASM)'
	@echo "Merging all per-sample VCFs..."
	@awk -F',' 'NR>1 && NF {gsub(/\r/,""); printf "results/%s/%s.vcf.gz\n",$$2,$$2}' $(DESIGN) > $(VCF_LIST)
	@$(MAKE) multisample
	@echo "All done. Merged VCF: $(MERGED_VCF)"

# -------------------------------------
# Merge VCFs
# -------------------------------------
vcf_list:
	@mkdir -p $(MERGED_DIR)
	@awk -F',' 'NR>1 && NF {gsub(/\r/,""); printf "results/%s/%s.vcf.gz\n",$$2,$$2}' $(DESIGN) > $(VCF_LIST)
	@echo "VCFs:"
	@cat $(VCF_LIST)

multisample: vcf_list
	@$(BCFTOOLS) merge -l $(VCF_LIST) -Oz -o $(MERGED_VCF)
	@$(TABIX) -p vcf $(MERGED_VCF)

# -------------------------------------
# Single sample helper
# -------------------------------------
one:
	@$(MAKE) all SRR=$(SRR) SAMPLE=$(SAMPLE) LAYOUT=$(LAYOUT) THREADS=$(THREADS) SUBSET=$(SUBSET) REF_ASM=$(REF_ASM)

# -------------------------------------
# Utility
# -------------------------------------
_ensure_dirs:
	@mkdir -p $(READS_DIR) $(GENOME_DIR) $(MERGED_DIR)

clean:
	@rm -rf $(READS_DIR) results $(GENOME_DIR)
	@rm -f $(GENOME).fai $(CHROMSIZES)