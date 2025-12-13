## Week 14: Perform an RNA-Seq differential gene expression study

### Overview

This project extends the Week 13 RNA-Seq workflow from read alignment and quantification to **differential expression analysis and functional enrichment**. Using the chr22 subset from the *Biostar Handbook*, I compared **Universal Human Reference (UHR)** RNA to **Human Brain Reference (HBR)** RNA to identify genes and biological processes that differ between the two reference transcriptomes.

The full analysis pipeline is automated with a Makefile and runs from a gene-level count matrix through visualization and enrichment.

---

## Input Data

* **Counts**: `counts/counts.csv` generated from featureCounts (Week 13)
* **Design file**: `design.csv` defining two groups (HBR vs UHR) with three replicates each

---

## Differential Expression Analysis

Differential expression was performed using **edgeR**, which models RNA-Seq count data using a negative binomial framework and accounts for biological variability across replicates.

The analysis produced a differential expression table (`results/edger.csv`) containing:

* log fold changes
* statistical significance (P-values)
* multiple-testing–corrected FDR values

Genes with **FDR < 0.05** were considered differentially expressed and carried forward for visualization and functional analysis.

---

## PCA Visualization

A PCA plot generated from the edgeR results shows **clear separation between HBR and UHR samples along PC1**, which explains **~97% of the variance**.

* Replicates cluster tightly within groups
* HBR and UHR samples are well separated
* This indicates strong biological signal and consistent sample processing

Overall, the PCA confirms that transcriptomic differences between UHR and HBR dominate the variance structure of the data.

---

## Heatmap of Differentially Expressed Genes

A heatmap of differentially expressed genes reveals two major gene clusters:

* Genes **up-regulated in UHR** and down-regulated in HBR
* Genes **up-regulated in HBR** and down-regulated in UHR

---

## Functional Enrichment Analysis

Genes with **FDR < 0.05** were subjected to Gene Ontology enrichment analysis using **g:Profiler**.

Functional enrichment of 288 differentially expressed genes identified 372 significant GO terms, dominated by metabolic, transport, and membrane-related processes.These results are consistent with known biological differences between a pooled cancer-cell RNA reference (UHR) and brain-derived RNA (HBR).

Enrichment results are provided in:

```
results/gprofiler.csv
```

---

## Reproducibility and Automation

The full analysis is automated using a Makefile with the following main steps:

* Convert featureCounts output to a count matrix
* Perform differential expression (edgeR)
* Generate PCA and heatmap visualizations
* Run functional enrichment analysis

Running:

```bash
make
```

reproduces all results from the count matrix onward.

---

## Summary

This project demonstrates a complete RNA-Seq analysis pipeline:

* robust differential expression between UHR and HBR
* strong agreement between statistical results and visualizations
* biologically interpretable functional enrichment

Together, these results validate both the computational workflow and the biological signal present in the data.