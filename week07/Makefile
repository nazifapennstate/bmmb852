# -------------------------------------
# Variables (overridable)
# -------------------------------------
SRR          ?= SRR1734993          # e.g., make all SRR=SRR10769653 LAYOUT=SE PLATFORM=ont
REF_ASM      ?= GCF_000848505.1
THREADS      ?= 4
SUBSET       ?= 939                  # set SUBSET= to fetch full dataset
LAYOUT       ?= PE                   # PE or SE
PLATFORM     ?= illumina             # illumina or ont

# Directories
GENOME_DIR   = genome
READS_DIR    = reads
RESULTS_DIR  = results/$(SRR)
SRA_CACHE    = $(HOME)/ncbi/public/sra

# -------------------------------------
# Files
# -------------------------------------
GENOME       = $(GENOME_DIR)/ebola_genome.fna
GFF          = $(GENOME_DIR)/ebola.gff
BAM          = $(RESULTS_DIR)/$(SRR)_sorted.bam
BAI          = $(BAM).bai
STATS        = $(RESULTS_DIR)/$(SRR)_alignment_stats.txt
BW           = $(RESULTS_DIR)/$(SRR)_coverage.bw
BG           = $(RESULTS_DIR)/$(SRR)_coverage.bedgraph
CHROMSIZES   = $(GENOME).chrom.sizes
MMI          = $(GENOME).mmi

# FASTQ naming for both layouts
READ1        = $(READS_DIR)/$(SRR)_1.fastq
READ2        = $(READS_DIR)/$(SRR)_2.fastq
READ_SE      = $(READS_DIR)/$(SRR).fastq

ifeq ($(LAYOUT),PE)
  FASTQ_TARGET = $(READ1) $(READ2)
else
  FASTQ_TARGET = $(READ_SE)
endif

# -------------------------------------
# Tools
# -------------------------------------
BWA          = bwa
MINIMAP2     = minimap2
SAMTOOLS     = samtools
BEDTOOLS     = bedtools
BG2BW        = bedGraphToBigWig
FASTQC       = fastqc
PREFETCH     = prefetch
FASTQ_DUMP   = fastq-dump
FASTERQ_DUMP = fasterq-dump
DATASETS     = datasets
UNZIP        = unzip

.PHONY: all fastq qc index mmi align stats bigwig clean

# -------------------------------------
# Default target
# -------------------------------------
all: qc stats bigwig

# -------------------------------------
# 1) Download genome (FASTA + optional GFF)
# -------------------------------------
$(GENOME):
	@echo "Downloading reference genome for $(REF_ASM)..."
	@mkdir -p $(GENOME_DIR)
	@$(DATASETS) download genome accession $(REF_ASM) \
	  --include genome,gff3 \
	  --filename $(GENOME_DIR)/$(REF_ASM).zip
	@$(UNZIP) -o $(GENOME_DIR)/$(REF_ASM).zip -d $(GENOME_DIR)/tmp_unzip
	@mv $(GENOME_DIR)/tmp_unzip/ncbi_dataset/data/*/*.fna $(GENOME)
	@if compgen -G "$(GENOME_DIR)/tmp_unzip/ncbi_dataset/data/*/*.gff" > /dev/null; then \
	  mv $(GENOME_DIR)/tmp_unzip/ncbi_dataset/data/*/*.gff $(GFF); \
	fi
	@rm -rf $(GENOME_DIR)/tmp_unzip $(GENOME_DIR)/$(REF_ASM).zip

# -------------------------------------
# 2) Download reads (PE or SE), subset optional
# -------------------------------------
fastq: $(FASTQ_TARGET)

$(FASTQ_TARGET):
	@echo "Downloading $(SRR) (LAYOUT=$(LAYOUT), SUBSET='$(SUBSET)')..."
	@mkdir -p $(READS_DIR)
	@$(PREFETCH) -O $(SRA_CACHE) $(SRR) || true
	@if [ -s "$(SRA_CACHE)/$(SRR).sra" ]; then SRC="$(SRA_CACHE)/$(SRR).sra"; else SRC="$(SRR)"; fi; \
	if [ "$(LAYOUT)" = "PE" ]; then \
	  if [ -n "$(SUBSET)" ]; then \
	    $(FASTQ_DUMP) -X $(SUBSET) --split-files -O $(READS_DIR) "$$SRC" || $(FASTQ_DUMP) --split-files -O $(READS_DIR) "$$SRC"; \
	  else \
	    $(FASTQ_DUMP) --split-files -O $(READS_DIR) "$$SRC" || $(FASTERQ_DUMP) --split-files -S -e $(THREADS) -O $(READS_DIR) "$$SRC"; \
	  fi; \
	else \
	  if [ -n "$(SUBSET)" ]; then \
	    $(FASTQ_DUMP) -X $(SUBSET) -O $(READS_DIR) "$$SRC" || $(FASTQ_DUMP) -O $(READS_DIR) "$$SRC"; \
	  else \
	    $(FASTQ_DUMP) -O $(READS_DIR) "$$SRC" || $(FASTERQ_DUMP) -S -e $(THREADS) -O $(READS_DIR) "$$SRC"; \
	  fi; \
	fi

# -------------------------------------
# 3) QC (works for both layouts)
# -------------------------------------
qc: fastq
	@echo "Running FastQC..."
	@mkdir -p $(READS_DIR)
	@inputs=""; \
	if [ -f "$(READ1)" ]; then inputs="$$inputs $(READ1)"; fi; \
	if [ -f "$(READ2)" ]; then inputs="$$inputs $(READ2)"; fi; \
	if [ -f "$(READ_SE)" ]; then inputs="$$inputs $(READ_SE)"; fi; \
	if [ -n "$$inputs" ]; then $(FASTQC) -t $(THREADS) -o $(READS_DIR) $$inputs; else echo "No FASTQ files found for $(SRR)"; exit 1; fi

# -------------------------------------
# 4) Index reference (BWA; harmless for ONT)
# -------------------------------------
index: $(GENOME)
	@echo "Indexing genome for BWA..."
	@$(BWA) index $(GENOME)

# Optional: minimap2 index for speed on ONT
mmi: $(GENOME)
	@echo "Building minimap2 index..."
	@$(MINIMAP2) -d $(MMI) $(GENOME)

# -------------------------------------
# 5) Align (auto-select by PLATFORM)
# -------------------------------------
align: index fastq
	@echo "Aligning $(SRR) (PLATFORM=$(PLATFORM), LAYOUT=$(LAYOUT))..."
	@mkdir -p $(RESULTS_DIR)
ifeq ($(PLATFORM),ont)
	@$(SAMTOOLS) faidx $(GENOME) >/dev/null 2>&1 || true
	@([ -f "$(MMI)" ] || $(MINIMAP2) -d $(MMI) $(GENOME))
	@$(MINIMAP2) -t $(THREADS) -ax map-ont $(MMI) $(READ_SE) \
	  | $(SAMTOOLS) sort -o $(BAM)
else
	@$(BWA) mem -t $(THREADS) $(GENOME) $(READ1) $(READ2) \
	  | $(SAMTOOLS) sort -o $(BAM)
endif
	@$(SAMTOOLS) index $(BAM)

# -------------------------------------
# 6) Alignment statistics
# -------------------------------------
stats: align
	@echo "Calculating alignment statistics for $(SRR)..."
	@mkdir -p $(RESULTS_DIR)
	@$(SAMTOOLS) flagstat $(BAM) > $(STATS)
	@$(SAMTOOLS) depth $(BAM) \
	  | awk '{s+=$$3; n++} END{if(n>0) printf "Average coverage = %.3f\n", s/n; else print "Average coverage = 0"}' \
	  >> $(STATS)
	@echo "Alignment stats written to $(STATS)"

# -------------------------------------
# 7) BigWig coverage
# -------------------------------------
bigwig: align
	@echo "Creating BigWig coverage track for $(SRR)..."
	@mkdir -p $(RESULTS_DIR)
	@$(SAMTOOLS) faidx $(GENOME)
	@cut -f1,2 $(GENOME).fai > $(CHROMSIZES)
	@$(BEDTOOLS) genomecov -ibam $(BAM) -bg \
	  | sort -k1,1 -k2,2n > $(BG)
	@$(BG2BW) $(BG) $(CHROMSIZES) $(BW)
	@echo "BigWig written to: $(BW)"

# -------------------------------------
# 8) Cleanup
# -------------------------------------
clean:
	@echo "Removing generated files and indices..."
	@rm -rf $(READS_DIR) $(RESULTS_DIR) $(GENOME_DIR)
	@rm -f $(GENOME).fai $(CHROMSIZES) $(MMI)
