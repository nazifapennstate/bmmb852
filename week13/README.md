## Week 13: Generate an RNA-Seq count matrix

For this assignment, I created an automated RNA-Seq workflow comparing **Universal Human Reference (UHR)** RNA vs **Human Brain Reference (HBR)** RNA using the chr22 subset from the *Biostar Handbook*. 

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
## 3. What the Makefile does

Running:

```bash
make
```

executes the full pipeline: alignment → BigWig coverage → gene-level count matrix.

| Command          | What It Does                        |
| ---------------- | ----------------------------------- |
| `make`           | Run **entire** pipeline end-to-end  |
| `make toolcheck` | Ensure required tools exist         |
| `make align`     | Produce BAM + BAI files             |
| `make bigwig`    | Produce .bw coverage tracks for IGV |
| `make counts`    | Generate the RNA-Seq count matrix   |
| `make clean`     | Remove output directories           |

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

<img width="841" height="184" alt="Screenshot 2025-12-09 at 11 20 28 PM" src="https://github.com/user-attachments/assets/d6956125-b0c4-45db-b8f3-709e76ad5e02" />

Many early chr22 genes show 0 counts in all samples, meaning these loci simply aren’t expressed in UHR or HBR. Farther down the matrix, some genes show higher counts in UHR than HBR, which matches what we see in the IGV coverage tracks.

## 5. IGV Visualisation

<img width="1431" height="436" alt="image" src="https://github.com/user-attachments/assets/90c55a0e-334b-4fc5-a501-b0fd67cd8fd7" />

* **Exon-restricted coverage patterns** → confirms RNA-Seq
* **Higher peaks in UHR vs HBR for some genes** → matches count matrix
* **Sharp intron dropouts** → typical of spliced mRNA

## 6. Using This Makefile With Other RNA-Seq Data

This Makefile can be reused for **any** single-end RNA-Seq dataset as long as:

1. **design.csv** is updated with new sample names
2. FASTQ files are placed in `reads/` and follow the naming pattern:

   ```
   <sample>_R1.fq
   ```
3. The **reference genome (REF)** and **annotation (GTF)** variables are changed to point to the correct organism or genomic region.
