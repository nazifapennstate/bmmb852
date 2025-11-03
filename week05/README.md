# Week 5: Obtain and visualise FASTQ data from SRA

I was assigned to Group 1, for which the paper was Gire et al. (2014), "Genomic surveillance elucidates Ebola virus origin and transmission during the 2014 outbreak," published in  [_Science_](https://doi.org/10.1126/science.1259657).

As outlined in the Acknowledgements section of the paper, the sequence data for this project are available at NCBI (NCBI BioGroup: PRJNA257197).

Under this bioproject, I picked a random strain (BioSample: SAMN03254248; Sample name: G5844.1; SRA: SRS803773) for all downstream analysis. 

## Writing a bash script

The bash script in scripts/download_illumina_ebola.sh does the following after activating the conda environment 'bioinfo':
1. Look up the sequencing runs for the SRA sample SRS803773.
2. Pick one run (paired-end if available).
3. Download a tiny test batch to measure the actual read length.
4. Calculate how many reads are needed to give ~10× coverage of the ~19 kb Ebola genome.
5. Download just that subset of reads (instead of the entire dataset).
6. Run quality checks: 
      - Count reads, total bases, and average read length (seqkit or awk).
      - Generate FastQC reports (HTML files with quality metrics).
7. Clean up the pilot test data and print a summary of what was downloaded.

### Coverage calculation in this script

In the Sequence Read Archive (SRA), a spot is the basic sequencing unit. For single-end runs, a spot contains one read, while for paired-end runs, a spot contains a read pair (R1 + R2). This matters because it determines how many bases each spot contributes.

To calculate how much sequencing data was needed, we started with the size of the Ebola genome, about 18,959 base pairs (~19 kb), and aimed for 10× coverage. That equals roughly 190,000 bases in total.

Each spot contributes a certain number of bases depending on read length. For paired-end runs, this means adding the average length of R1 and R2. In this dataset, the pilot download showed that each read was about 101 base pairs long, so each spot contributed about 202 bases.

Dividing the total number of bases needed by the amount contributed per spot gave 939 spots, which the script then downloaded. These appear as 939 reads in R1 and 939 reads in R2, providing the target 10× coverage of the Ebola genome.

## Quality assessment

Generate basic statistics on the downloaded reads (e.g., number of reads, total bases, average read length).

```bash
file                      format  type  num_seqs  sum_len  min_len  avg_len  max_len
reads/SRR1734993_1.fastq  FASTQ   DNA        939   94,839      101      101      101
reads/SRR1734993_2.fastq  FASTQ   DNA        939   94,839      101      101      101
```

Run FastQC on the downloaded data to generate a quality report.

<img width="1438" height="803" alt="image" src="https://github.com/user-attachments/assets/0a63c177-4d31-4fe0-bff1-f3bf462f4630" />

The FastQC report shows quality drops at the read ends, suggesting low-quality bases could be trimmed.

## Compare sequencing platforms 

I searched for Ebola virus genomes sequenced by Nanopore and chose a random strain (SRR10769653). The script for this can be found in scripts/download_nanopore_ebola.sh.

```bash
file                     format  type  num_seqs  sum_len  min_len  avg_len  max_len
reads/SRR10769653.fastq  FASTQ   DNA        300   77,650       35    258.8      695
```

<img width="1434" height="796" alt="image" src="https://github.com/user-attachments/assets/368236aa-f562-4720-81dd-61e4b967ff8d" />

Compared to the Illumina dataset, the Nanopore dataset shows much lower per-base quality.
