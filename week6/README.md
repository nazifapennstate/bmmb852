# Week 6: Generate a BAM alignment file

This week’s goal was to transform the Week 5 Bash pipeline into a reproducible Makefile which does the following:

| Target                                     | Purpose                                                                         |
| ------------------------------------------ | ------------------------------------------------------------------------------- |
| `all`                                      | Default target → runs the full pipeline through `stats`                         |
| `$(GENOME)`                                | Downloads and extracts the Ebola reference genome (RefSeq GCF_000848505.1)      |
| `$(READS_DIR)/$(SRR)_1.fastq` + `_2.fastq` | Downloads a 10× subset of reads from SRA (SRR1734993)                           |
| `qc`                                       | Runs FastQC on downloaded reads                                                 |
| `index`                                    | Builds the BWA index for the reference genome                                   |
| `align`                                    | Aligns reads → generates sorted BAM + index (`samtools sort`, `samtools index`) |
| `stats`                                    | Summarizes alignment statistics and mean coverage                               |
| `clean`                                    | Removes all intermediate directories and temporary files                        |

After activating the bioinfo Conda environment, the entire analysis can be reproduced with one command:

```bash
make all
```
Individual steps (e.g., make qc or make align) can also be run independently. All resulting files are stored in structured subdirectories:

```bash
genome/   → reference FASTA + index files
reads/    → paired FASTQ files + FastQC reports
results/  → BAM + BAI + alignment statistics
```
## Visualising BAM file

<img width="1305" height="405" alt="Screenshot 2025-10-05 at 8 02 58 PM" src="https://github.com/user-attachments/assets/56ab0573-1b13-488d-83a7-040b91366cf0" />

## Alignment statistics

From the output results/alignment_stats.txt:
```bash
1881 + 0 in total (QC-passed reads + QC-failed reads)
1878 + 0 primary
32 + 0 mapped (1.70%)
29 + 0 primary mapped (1.54%)
26 + 0 properly paired (1.38%)
Average coverage = 1.61×
```
### _What percentage of reads aligned to the genome?_
Approximately 1.7% of reads aligned to the genome. 

### _What was the expected average coverage?_
The expected coverage was ~10×, based on downloading a subset (~939 paired-end reads) to achieve roughly tenfold coverage of the ~19 kb Ebola genome.

### _What is the observed average coverage?_
The observed coverage was 1.61×, calculated using samtools depth in the Makefile.

### _How much does the coverage vary across the genome?_
Coverage varied dramatically across the genome. In IGV, only small regions showed mapped reads, while most of the genome remained uncovered. This unevenness likely reflects that the 10× subset we downloaded achieved only ~1.7× actual coverage, leaving many genomic regions with too few reads to align.
