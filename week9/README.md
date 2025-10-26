# Week 9: Revising and Improving Your Automation Code

This week I took the automated workflow I built in Week 8 and made it easier to run and debug.
I wanted something I could execute in **one command** (`make run`) that checks tools, fetches the reference if missing, and launches a parallelized batch run over a small subset of samples for quick testing.

## 1. Building the `design.csv` file

I started again from BioProject **PRJNA257197** (Gire et al., 2014 *Science*, Ebola outbreak).
I downloaded the full run metadata and then filtered it down to **one SRR per unique sample**, keeping only **five samples** for speed.

```bash
# 1. Fetch SRA RunInfo
esearch -db sra -query PRJNA257197 | efetch -format runinfo > sra_runinfo.csv

# 2. Extract one SRR per unique Sample (first 5)
awk -F, 'BEGIN{OFS=","}
NR==1{
  for(i=1;i<=NF;i++) h[$i]=i
  print "Run,Sample,Layout"; next
}
{
  run=$h["Run"];
  sample=$h["SampleName"]; if(sample=="") sample=$h["Sample"]; if(sample=="") sample=$h["BioSample"];
  layout=toupper($h["LibraryLayout"]); if(layout ~ /PAIRED|PE/) layout="PE"; else layout="SE";
  if(run!="" && sample!="" && !seen[sample]++) print run,sample,layout
}' sra_runinfo.csv | head -n 6 > design.csv
```

Each line of `design.csv` contains:

* **Run** → the SRA accession (SRR)
* **Sample** → the library or patient ID
* **Layout** → `PE` or `SE` depending on library type

---

## 2. Updating the Makefile

I reused my Week 8 Makefile but streamlined it:

* Deleted the inline comments on variable definitions.
* Added a `run` target that performs:

  1. a dependency check for required tools,
  2. genome setup if missing,
  3. and then executes the batch analysis using GNU Parallel.
* Kept the detailed `echo` statements from earlier weeks so I can see exactly what’s happening at each step—FastQC, BWA MEM alignment, BAM sorting, coverage, and BigWig generation.

The structure still uses the same logical targets:
`setup-genome → fastq → qc → align → stats → bigwig`.

---

## 3. Running the workflow

To test everything on my five-sample design file, I used:

```bash
make run DESIGN=design.csv JOBS=2 THREADS=2 SUBSET=10000
```

Here’s what happens when I run this command:

* It checks that all required tools (`bwa`, `samtools`, `bedtools`, `fastqc`, etc.) are installed.
* If the genome isn’t present, it downloads and indexes it automatically.
* Then it launches **two jobs in parallel**, each using two threads.
* `SUBSET=10000` limits the analysis to the first 10,000 reads of each SRR—enough to validate the pipeline quickly.

---

## 4. Comparing `make run`, `make batch`, and `make all`

* **`make run`** — It wraps the entire workflow: checking tools, verifying genome files, and then calling `make batch`.
  It’s the easiest way to execute the full pipeline from scratch.

* **`make batch`** — Runs the analysis for **all rows** in `design.csv` using GNU Parallel.
  Each row’s variables (`SRR`, `Sample`, `Layout`) are passed to `make all`.

* **`make all`** — Runs the workflow for **a single sample** (FastQC → alignment → stats → BigWig).
  It can be used to troubleshoot a single SRR if something goes wrong in the batch.

In summary, `run` handles the setup, `batch` loops through samples, and `all` runs one sample.

---

## 5. Output structure

Each sample’s results are organized in the same structure as before:

```bash
genome/                       # reference FASTA + index
reads/                        # FASTQs + FastQC reports
results/<Sample>/
  ├── <Sample>_sorted.bam
  ├── <Sample>_sorted.bam.bai
  ├── <Sample>_alignment_stats.txt
  ├── <Sample>_coverage.bedgraph
  └── <Sample>_coverage.bw
logs/                         # stdout/stderr for each SRR
parallel.log                  # summary table of job status
```

---

## 6. Tools used

* **bwa** — sequence alignment
* **samtools** — sorting, indexing, and coverage depth
* **bedtools** — convert BAM depth to bedGraph
* **bedGraphToBigWig** — generate browser tracks
* **fastqc** — read quality control
* **parallel** — batch execution
* **datasets / entrez-direct** — genome and metadata retrieval
* **sra-tools** — download raw reads (`prefetch`, `fastq-dump`, `fasterq-dump`)

---

In short, this week I focused on **cleaning up my Makefile**, **removing unnecessary inline comments**, and creating a single-command, reproducible workflow.
The `make run` command now serves as a convenient entry point for large-scale or quick-subset analyses alike.