# ============================================================
# Makefile: RNA-Seq DE + enrichment (builds on Week 13 outputs)
# ============================================================

SHELL := /usr/bin/env bash
.SHELLFLAGS := -euo pipefail -c

THREADS ?= 4
DESIGN  ?= design.csv

COUNT_DIR := counts
RES_DIR   := results

FC_TXT := $(COUNT_DIR)/counts.txt
COUNTS := $(COUNT_DIR)/counts.csv

# Choose edger.csv or deseq2.csv as your canonical result
DE_RES := $(RES_DIR)/edger.csv
# DE_RES := $(RES_DIR)/deseq2.csv

PCA_PDF     := $(RES_DIR)/pca.pdf
HEATMAP_PDF := $(RES_DIR)/heatmap.pdf

GP_CSV := $(RES_DIR)/gprofiler.csv

.PHONY: all dirs counts_csv de pca heatmap gprofiler clean

all: $(DE_RES) $(PCA_PDF) $(HEATMAP_PDF) $(GP_CSV)

dirs:
	mkdir -p $(COUNT_DIR) $(RES_DIR)

# --- Convert featureCounts output to counts.csv (stats env) ---
$(COUNTS): $(FC_TXT) | dirs
	micromamba run -n stats Rscript src/r/format_featurecounts.r -c $(FC_TXT) -o $@

counts_csv: $(COUNTS)

# --- Differential expression (stats env) ---
$(RES_DIR)/edger.csv: $(COUNTS) $(DESIGN) | dirs
	micromamba run -n stats Rscript src/r/edger.r -d $(DESIGN) -c $(COUNTS)
	mv -f edger.csv $@

$(RES_DIR)/deseq2.csv: $(COUNTS) $(DESIGN) | dirs
	micromamba run -n stats Rscript src/r/deseq2.r -d $(DESIGN) -c $(COUNTS)
	mv -f deseq2.csv $@

de: $(DE_RES)

# --- PCA + heatmap (stats env) ---
$(PCA_PDF): $(DE_RES) $(DESIGN) | dirs
	micromamba run -n stats Rscript src/r/plot_pca.r -c $(DE_RES) -d $(DESIGN) -o $@

$(HEATMAP_PDF): $(DE_RES) $(DESIGN) | dirs
	micromamba run -n stats Rscript src/r/plot_heatmap.r -c $(DE_RES) -d $(DESIGN) -o $@

pca: $(PCA_PDF)
heatmap: $(HEATMAP_PDF)

# --- Enrichment (bioinfo env) ---
$(RES_DIR)/gprofiler.csv: $(DE_RES) | dirs
	micromamba run -n bioinfo bio gprofiler -c $(DE_RES) -d hsapiens
	mv -f gprofiler.csv $@

gprofiler: $(GP_CSV)

clean:
	rm -rf $(RES_DIR) $(COUNTS)