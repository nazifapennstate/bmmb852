# Week 8: Automate a large scale analysis

This workflow uses sequencing data from the Gire et al. (2014) *Science* paper on the **2014 West African Ebola outbreak** (BioProject **PRJNA257197**). 

## 1. Building the design.csv file

To start, I generated `sra_runinfo.csv` for BioProject **PRJNA257197** using NCBI’s Entrez Direct utilities (`esearch`, `efetch`):

```bash
esearch -db sra -query PRJNA257197 | efetch -format runinfo > sra_runinfo.csv
````

Then I used `awk` to extract and clean up the key columns into a consistent format:

```bash
awk -F, 'BEGIN{OFS=","}
NR==1{
  for(i=1;i<=NF;i++) h[$i]=i
  print "Run,Sample,Layout"; next
}
{
  run=$h["Run"];
  sample=$h["SampleName"]; if(sample=="") sample=$h["Sample"]; if(sample=="") sample=$h["BioSample"];
  layout=toupper($h["LibraryLayout"]); if(layout!="PE" && layout!="SE") layout="PE";
  if(run!="" && sample!="") print run,sample,layout
}' sra_runinfo.csv | LC_ALL=C sort -u > design.csv
```

### What each column means

* **Run** → SRA run accession (SRR) from the RunInfo table; each line represents one sequencing run.
* **Sample** → the submitter’s sample or library name, taken from `SampleName`, or from `Sample`/`BioSample` if missing.
* **Layout** → `LibraryLayout` (paired-end `PE` or single-end `SE`).

---

## 2. Reasoning about sample counts

When I counted the unique sample names in my `design.csv`, I found **713** unique IDs:

```bash
awk -F, 'NR>1 {print $2}' design.csv | sort -u | wc -l
# 713
```

At first this seems inconsistent with the paper, which reports sequencing **99 genomes from 78 patients**, but it appears that the SRA BioProject **PRJNA257197** contains all raw sequencing runs (713 unique libraries, 891 total SRRs), not just the final assemblies. This might include:

* Technical replicates.
* Repeated library preps from the same patient.
* Negative controls and test libraries.
* Failed or partial runs.

The paper states that 99 complete viral genomes were generated from 78 patients, but it does not list which specific SRRs were used. I could not find any metadata connecting individual SRRs to the final 99 genomes. Without that, I cannot know which replicate was used in the analysis.

For this assignment, I treat the full `design.csv` as the input list and let the Makefile process each row as its own sample.

---

## 3. Automating batch analysis

Once `design.csv` is ready, the Makefile can be used to run the full workflow. The Makefile and design.csv should be placed in the same working directory before running the commands. 

_Note on Makefile modification: This Makefile is adapted from my Week 7 workflow, where the Makefile supported both Illumina (BWA MEM) and Oxford Nanopore (minimap2) reads. Since the Gire et al. (2014) dataset contains only Illumina reads, I simplified the Week 7 version by removing the minimap2 steps and related variables. The rest of the structure—indexing, alignment, statistics, and BigWig generation—remains consistent with Week 7 for continuity and reproducibility._

### Step 1. Download and index reference genome

```bash
make setup-genome REF_ASM=GCF_000848505.1 THREADS=4
```

### Step 2. Preview commands (dry run)

```bash
make batch DESIGN=design.csv JOBS=4 THREADS=4 DRYRUN=1
```

### Step 3. Run the full workflow

```bash
make batch DESIGN=design.csv JOBS=4 THREADS=4
```

Optional quick testing with a smaller subset:

```bash
make batch DESIGN=design.csv JOBS=4 THREADS=4 SUBSET=100000
```
_Note: setting SUBSET=100000 means that only the first 100,000 reads are retrieved from each SRA run_

---

## 4. Output structure

Each sample produces its own directory under `results/`:

```bash
genome/                       # reference FASTA + index
reads/                        # sample FASTQs + FastQC reports
results/<Sample>/
  ├── <Sample>_sorted.bam
  ├── <Sample>_sorted.bam.bai
  ├── <Sample>_alignment_stats.txt
  ├── <Sample>_coverage.bedgraph
  └── <Sample>_coverage.bw
```

---

## 5. Tools used

* **bwa** — sequence alignment
* **samtools** — sorting, indexing, coverage depth
* **bedtools** — coverage to bedGraph
* **bedGraphToBigWig** — convert coverage to browser track
* **fastqc** — read quality control
* **parallel** — concurrent sample execution
* **datasets / entrez-direct** — genome and metadata retrieval
* **sra-tools** — downloading raw reads (`prefetch`, `fastq-dump`, `fasterq-dump`)
