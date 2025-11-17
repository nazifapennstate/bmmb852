# Week 12: Evaluate data from the Cancer Genome in a Bottle project

For this assignment, I chose option 1 and evaluated alignments across sequencing platforms. I automated everything with a Makefile. 

Running:

```bash
make all
````

executes the full workflow.

---

## What the Makefile does

Running `make all` performs:

1. **Download of the GRCh38 GIABv3 reference genome** and indexing with `samtools faidx`.

2. **Direct HTTPS extraction of BAM alignments** for:

   * Illumina WGS (NYGC, 118×)
   * Element AVITI (77×)
   * PacBio HiFi Revio (35×)

3. **Subsetting to the KRAS gene region**:

   ```
   chr12:25205246–25250936
   ```

4. **Creation of per-platform BAM files**, indexing, and generation of:

   * `*.flagstat.txt` (mapping statistics)
   * `*.depth.txt` (coverage at each base)

Outputs are written to the `bam/` directory.

---

## How statistics were calculated

### 1. Mapping statistics

Using `samtools flagstat`:

```bash
samtools flagstat bam/kras_illumina.bam > bam/kras_illumina.flagstat.txt
```

Fields interpreted:

* total reads in region
* % mapped
* % properly paired (short-read platforms)

### 2. Coverage

Mean depth across the region was calculated with:

```bash
awk '{sum+=$3} END {print sum/NR}' bam/kras_platform.depth.txt
```

This computes the average of the depth column for the entire KRAS window.

---

## Summary of Results

### Coverage (Average Depth)

| Platform          | Avg. Depth (×) |
| ----------------- | -------------- |
| **Illumina**      | **118.9×**     |
| **Element AVITI** | **89.9×**      |
| **PacBio HiFi**   | **45.6×**      |

### Mapping / Pairing (from `flagstat`)

| Platform          | Total Reads | % Mapped | Properly Paired |
| ----------------- | ----------- | -------- | --------------- |
| **Illumina**      | 45,834      | 99.78%   | 98.43%          |
| **Element AVITI** | 28,035      | 99.85%   | 99.25%          |
| **PacBio HiFi**   | 149         | 100%     | N/A             |

---

## Alignment visualisation in IGV

Below are representative IGV screenshots from the KRAS region
(`chr12:25205246–25250936`) for each sequencing platform:

### **Illumina (NYGC 118×)**
<img width="2048" height="1002" alt="image" src="https://github.com/user-attachments/assets/6ad47090-9d7b-41c1-867e-19ebc23d309d" />

### **Element AVITI (77×)**
<img width="2048" height="1210" alt="image" src="https://github.com/user-attachments/assets/f55fcf8b-cc37-4593-8c95-6fe620986765" />

### **PacBio HiFi Revio (35×)**
<img width="2048" height="973" alt="image" src="https://github.com/user-attachments/assets/965c20d7-cc8c-488f-b9d8-cc9c99063820" />

---

## Interpretation

Illumina provided the highest depth, Element showed the strongest short-read alignment metrics, and PacBio HiFi produced the most accurate long-read alignments despite lower coverage.
