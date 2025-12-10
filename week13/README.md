## Week 13: Generate an RNA-Seq count matrix

For this assignment, I created an automated RNA-Seq workflow comparing **Universal Human Reference (UHR)** RNA vs **Human Brain Reference (HBR)** RNA using the chr22 subset from the *Biostar Handbook*. Running:

```bash
make
```

executes the full pipeline: alignment → BigWig coverage → gene-level count matrix.

---

## 1. Downloading the Data

The dataset is a curated subset containing only **chromosome 22** plus annotations. Download and unpack with:

```bash
wget -nc http://data.biostarhandbook.com/data/uhr-hbr.tar.gz
tar xzvf uhr-hbr.tar.gz
```

This creates:

```
refs/chr22.genome.fa
refs/chr22.gtf
reads/HBR_*.fq
reads/UHR_*.fq
```

### What the samples are

* **UHR (Universal Human Reference)** → total RNA from *10 cancer cell lines*
* **HBR (Human Brain Reference)** → total RNA from *23 human brain donors*

Both are commercially produced reference RNAs used widely for benchmarking RNA-Seq workflows.

---

## 2. Design File

Create a file named `design.csv`:

```
sample,group
HBR_1,HBR
HBR_2,HBR
HBR_3,HBR
UHR_1,UHR
UHR_2,UHR
UHR_3,UHR
```

The Makefile uses the `sample` column to locate FASTQ files and generate outputs automatically.

---

## 3. Pipeline Overview (Makefile)

### **Alignment**

HISAT2 aligns each sample to chr22, producing:

```
bam/<sample>.bam
bam/<sample>.bam.bai
```

### **Coverage Tracks**

BigWig files are generated for visualization in IGV:

```
bw/<sample>.bw
```

### **Gene-Level Counts**

FeatureCounts summarizes expression across all six samples:

```
counts/counts.txt
```

---

## 4. Understanding the Count Matrix

FeatureCounts produces columns:

| Column                      | Meaning                                        |
| --------------------------- | ---------------------------------------------- |
| **Geneid**                  | Ensembl gene ID                                |
| **Chr, Start, End, Strand** | Genomic gene coordinates                       |
| **Length**                  | Total exon length                              |
| **HBR / UHR columns**       | Number of reads overlapping exons of that gene |


```

→ gene expressed in UHR_1 and UHR_2, absent in HBR samples.

---

## 5. IGV Verification


* **Exon-restricted coverage patterns** → confirms RNA-Seq
* **Higher peaks in UHR vs HBR for some genes** → matches count matrix
* **Sharp intron dropouts** → typical of spliced mRNA