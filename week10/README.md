# Week 10: Generate a multisample variant call file (VCF)

This version of the Makefile:

* Accepts a `design.csv` file containing Run, Sample, and Layout columns.  
* Automatically downloads and indexes the reference genome if missing.  
* Runs **FastQC → BWA MEM → BAM sort → stats → BigWig → VCF** for each sample.  
* Produces a **merged multi-sample VCF** (`results/merged/merged.vcf.gz`).  
* Organizes all intermediate and final files under structured folders.  

Everything can be executed in **one command**:

```bash
make run DESIGN=design.csv JOBS=4 THREADS=4 SUBSET=10000
````

---

## 2. Why There Are Both `SRR` and `SAMPLE` Columns

The pipeline separates **technical runs** from **biological samples**:

| Variable     | Meaning                                                               | Example      |
| ------------ | --------------------------------------------------------------------- | ------------ |
| **`SRR`**    | SRA Run accession: represents an individual sequencing run or lane.  | `SRR1553418` |
| **`SAMPLE`** | Biological or clinical sample identifier: can include multiple SRRs. | `EM096`      |

This distinction ensures that if one sample has several sequencing runs (technical replicates), their data can later be grouped or merged under a single biological name.

For example:

| Run        | Sample | Layout |
| ---------- | ------ | ------ |
| SRR1553418 | EM096  | PE     |
| SRR1553419 | EM096  | PE     |

Here both SRRs map to the same biological sample (`EM096`).
During batch execution, results are stored in `results/EM096/`, preventing duplication or overwrites.

When testing a single SRR manually, the user can still assign the same name for both:

```bash
make all SRR=SRR1553418 SAMPLE=EM096 LAYOUT=PE THREADS=4 SUBSET=10000
```

If the user omits `SAMPLE=`, it defaults to the SRR name (e.g., `results/SRR1553418/`).

---

## 3. Debugging Coverage Issues

Initially, I ran the workflow on a few SRRs from the *Gire et al.*, 2014 *Science* Ebola dataset.
Some SRRs produced very low mapping percentages.

**Example of poor coverage (initial SRRs):**


After inspecting `*_alignment_stats.txt` and `samtools flagstat` results, I swapped out those low-coverage SRRs for better ones.

**Improved coverage (final dataset):**



## 4. Running the Pipeline

### Single sample

```bash
make all SRR=SRR1553418 SAMPLE=EM096 LAYOUT=PE THREADS=4 SUBSET=10000
```

### Batch mode (all samples in design.csv)

```bash
make run DESIGN=design.csv JOBS=4 THREADS=4 SUBSET=10000
```

### Clean up

```bash
make clean
```

---

## 5. Output Structure

```bash
genome/                     # reference FASTA + GFF + index
reads/                      # raw FASTQs + FastQC reports
results/<Sample>/           # per-sample results
  ├── <Sample>_sorted.bam
  ├── <Sample>_alignment_stats.txt
  ├── <Sample>_coverage.bw
  └── <Sample>.vcf.gz
results/merged/             # merged multi-sample VCF
  └── merged.vcf.gz
.tmp/                       # cleaned design.csv (temporary)
```

---

## 6. Tools Used

| Tool                   | Purpose                                                  |
| ---------------------- | -------------------------------------------------------- |
| **bwa**                | sequence alignment                                       |
| **samtools**           | BAM sorting, indexing, coverage                          |
| **bedtools**           | genome coverage → bedGraph                               |
| **bedGraphToBigWig**   | browser track generation                                 |
| **bcftools**           | variant calling + merging                                |
| **fastqc**             | read quality control                                     |
| **parallel**           | concurrent batch execution                               |
| **datasets / esearch** | genome retrieval                                         |
| **sra-tools**          | read download (`prefetch`, `fastq-dump`, `fasterq-dump`) |

---

## 7. Key Commands Recap

| Command             | Description                      |
| ------------------- | -------------------------------- |
| `make toolcheck`    | verify all tools installed       |
| `make setup-genome` | download + index reference       |
| `make all`          | run full workflow for one sample |
| `make batch`        | run all samples in design.csv    |
| `make multisample`  | merge per-sample VCFs            |
| `make clean`        | remove reads/, results/, genome/ |