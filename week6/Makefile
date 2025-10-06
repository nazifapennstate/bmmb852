# Set the shell the commands run in.
SHELL = bash

# Execute all commands in a single shell.
.ONESHELL:

# Run the shell with strict error checking.
.SHELLFLAGS = -eu -o pipefail -c

# Delete target files if the command fails.
.DELETE_ON_ERROR:

# Warn if a variable is not defined.
MAKEFLAGS += --warn-undefined-variables

# Disable built-in rules.
MAKEFLAGS += --no-builtin-rules


# -------------------------------------
# Variables
# -------------------------------------
# Ebola (Gire et al. 2014; RefSeq assembly for NC_002549.1)
SRS          = SRS803773
SRR          = SRR1734993
REF_ASM      = GCF_000848505.1
THREADS      = 4

# Organized directories
GENOME_DIR   = genome
READS_DIR    = reads
RESULTS_DIR  = results
SRA_CACHE    = $(HOME)/ncbi/public/sra

# ------ No edits beyond this point ------

# Files
GENOME       = $(GENOME_DIR)/ebola_genome.fna
BAM          = $(RESULTS_DIR)/ebola_sorted.bam
BAI          = $(RESULTS_DIR)/ebola_sorted.bam.bai
STATS        = $(RESULTS_DIR)/alignment_stats.txt

# Tools
BWA          = bwa
SAMTOOLS     = samtools
FASTQC       = fastqc
CURL         = curl
PREFETCH     = prefetch
FASTQ_DUMP   = fastq-dump
FASTERQ_DUMP = fasterq-dump
DATASETS     = datasets
UNZIP        = unzip

.PHONY: all qc index align stats clean

# -------------------------------------
# Default target
# -------------------------------------
all: stats

# -------------------------------------
# 1) Download genome (RefSeq assembly -> FASTA)
# -------------------------------------
$(GENOME):
	mkdir -p $(GENOME_DIR)
	$(DATASETS) download genome accession $(REF_ASM) \
	  --include genome,gff3 \
	  --filename $(GENOME_DIR)/$(REF_ASM).zip
	$(UNZIP) -o $(GENOME_DIR)/$(REF_ASM).zip -d $(GENOME_DIR)/tmp_unzip
	# Move FASTA to canonical path; keep GFF if present
	mv $(GENOME_DIR)/tmp_unzip/ncbi_dataset/data/*/*.fna $(GENOME)
	if compgen -G "$(GENOME_DIR)/tmp_unzip/ncbi_dataset/data/*/*.gff" > /dev/null; then \
	  mv $(GENOME_DIR)/tmp_unzip/ncbi_dataset/data/*/*.gff $(GENOME_DIR)/ebola.gff; \
	fi
	rm -rf $(GENOME_DIR)/tmp_unzip $(GENOME_DIR)/$(REF_ASM).zip

# -------------------------------------
# 2) Get SRA locally and make a 10× subset of reads
#    Primary: prefetch + fastq-dump -X 939 (paired)
#    Fallbacks: fastq-dump full; then fasterq-dump full
# -------------------------------------
$(READS_DIR)/$(SRR)_1.fastq $(READS_DIR)/$(SRR)_2.fastq:
	mkdir -p $(READS_DIR)
	# Ensure local cache SRA object exists
	$(PREFETCH) -O $(SRA_CACHE) $(SRR) || true
	if [ -s "$(SRA_CACHE)/$(SRR).sra" ]; then \
	  echo "Using local $(SRA_CACHE)/$(SRR).sra"; \
	  if $(FASTQ_DUMP) -X 939 --split-files -O $(READS_DIR) "$(SRA_CACHE)/$(SRR).sra"; then \
	    echo "Downloaded 10× subset (939 spots) via fastq-dump."; \
	  else \
	    echo "Subset failed; trying full fastq-dump..."; \
	    $(FASTQ_DUMP) --split-files -O $(READS_DIR) "$(SRA_CACHE)/$(SRR).sra"; \
	  fi; \
	else \
	  echo "No local .sra; trying direct fastq-dump (subset) from SRA..."; \
	  if $(FASTQ_DUMP) -X 939 --split-files -O $(READS_DIR) $(SRR); then \
	    echo "Downloaded 10× subset (939 spots) via fastq-dump (remote)."; \
	  else \
	    echo "Subset via fastq-dump failed; trying full fastq-dump..."; \
	    if $(FASTQ_DUMP) --split-files -O $(READS_DIR) $(SRR); then \
	      echo "Downloaded full dataset via fastq-dump."; \
	    else \
	      echo "fastq-dump failed; trying fasterq-dump (full) ..."; \
	      $(FASTERQ_DUMP) --split-files -S -e $(THREADS) -O $(READS_DIR) $(SRR); \
	    fi; \
	  fi; \
	fi

# -------------------------------------
# 3) Run FastQC
# -------------------------------------
qc: $(READS_DIR)/$(SRR)_1.fastq $(READS_DIR)/$(SRR)_2.fastq
	$(FASTQC) -o $(READS_DIR) $(READS_DIR)/$(SRR)_1.fastq $(READS_DIR)/$(SRR)_2.fastq

# -------------------------------------
# 4) Index the genome (BWA)
# -------------------------------------
index: $(GENOME)
	$(BWA) index $(GENOME)

# -------------------------------------
# 5) Align reads → sorted, indexed BAM
# -------------------------------------
align: index $(READS_DIR)/$(SRR)_1.fastq $(READS_DIR)/$(SRR)_2.fastq
	mkdir -p $(RESULTS_DIR)
	$(BWA) mem -t $(THREADS) $(GENOME) \
	  $(READS_DIR)/$(SRR)_1.fastq $(READS_DIR)/$(SRR)_2.fastq \
	| $(SAMTOOLS) sort -o $(BAM)
	$(SAMTOOLS) index $(BAM)

# -------------------------------------
# 6) Alignment statistics (flagstat + mean depth)
# -------------------------------------
stats: align
	$(SAMTOOLS) flagstat $(BAM) > $(STATS)
	$(SAMTOOLS) depth $(BAM) \
	  | awk '{s+=$$3; n++} END{if(n>0) printf "Average coverage = %.3f\n", s/n; else print "Average coverage = 0"}' \
	  >> $(STATS)
	@echo "Alignment stats written to $(STATS)"

# -------------------------------------
# 7) Cleanup
# -------------------------------------
clean:
	rm -rf $(GENOME_DIR) $(READS_DIR) $(RESULTS_DIR) *.csv *.html *.zip *.txt
