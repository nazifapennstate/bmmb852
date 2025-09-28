# Week 5: Obtain and visualise FASTQ data from SRA

I was assigned to Group 1, for which the paper was Gire et al. (2014), "Genomic surveillance elucidates Ebola virus origin and transmission during the 2014 outbreak," published in  [_Science_](https://doi.org/10.1126/science.1259657).

As outlined in the Acknowledgements section of the paper, the sequence data for this project are available at NCBI (NCBI BioGroup: PRJNA257197).

Under this bioproject, I picked a random strain (BioSample: SAMN03254248; Sample name: G5844.1; SRA: SRS803773) for all downstream analysis. 

## Writing a bash script

Below is a bash script that does the following after activating the conda environment 'bioinfo':
1. Look up the sequencing runs for the SRA sample SRS803773.
2. Pick one run (paired-end if available).
3. Download a tiny test batch to measure the actual read length.
4. Calculate how many reads are needed to give ~10× coverage of the ~19 kb Ebola genome.
5. Download just that subset of reads (instead of the entire dataset).
6. Run quality checks: 
      - Count reads, total bases, and average read length (seqkit or awk).
      - Generate FastQC reports (HTML files with quality metrics).
7. Clean up the pilot test data and print a summary of what was downloaded.

```bash
# Set flags for ease of coding
set -xeuo pipefail

# Set inputs
SRS=SRS803773          # SRA Sample
genome_size=18959      # bp (Ebola ~18,959)
coverage=10            # target depth

# ------ No edits beyond this point ------

# 1) Get SRA RunInfo for the sample
esearch -db sra -query "$SRS" | efetch -format runinfo > ${SRS}_runinfo.csv

# 2) Choose a run (prefer paired-end if available)
#    Run=col1, avgLength=col7, LibraryLayout=col16
SRR=$(awk -F',' 'NR==1{next} $16=="PAIRED"{print $1; exit}' ${SRS}_runinfo.csv)
[ -z "$SRR" ] && SRR=$(awk -F',' 'NR==1{next}{print $1; exit}' ${SRS}_runinfo.csv)

layout=$(awk -F',' -v r="$SRR" 'NR==1{next} $1==r {print $16; exit}' ${SRS}_runinfo.csv)

# 3) Pilot download: grab a small subset to measure REAL read length
mkdir -p pilot reads
fastq-dump -X 100 --split-files --outdir pilot "$SRR"

# 4) Measure bases per SPOT from the pilot 
if [ "$layout" = "PAIRED" ]; then
  # mean length of R1 and R2; bases per spot = mean(R1) + mean(R2)
  mean_r1=$(awk 'NR%4==2{n++; s+=length($0)} END{if(n>0) printf "%.0f", s/n; else print 0}' pilot/${SRR}_1.fastq)
  mean_r2=$(awk 'NR%4==2{n++; s+=length($0)} END{if(n>0) printf "%.0f", s/n; else print 0}' pilot/${SRR}_2.fastq)
  bases_per_spot=$(( mean_r1 + mean_r2 ))
else
  mean_r1=$(awk 'NR%4==2{n++; s+=length($0)} END{if(n>0) printf "%.0f", s/n; else print 0}' pilot/${SRR}.fastq)
  bases_per_spot=$(( mean_r1 ))
fi

# 5) Compute spots needed for ~10 coverage (ceil division)
target_bases=$(( coverage * genome_size ))
spots_needed=$(( (target_bases + bases_per_spot - 1) / bases_per_spot ))

# 6) Final download: only what you need
fastq-dump -X "$spots_needed" --split-files --outdir reads "$SRR"

# 7) QC: stats + FastQC
if command -v seqkit >/dev/null 2>&1; then
  seqkit stats reads/${SRR}*.fastq > reads/${SRR}_stats.txt
else
  # awk fallback if seqkit isnt installed
  {
    echo -e "file\treads\ttotal_bases\tavg_len"
    for f in reads/${SRR}*.fastq; do
      awk -v fn="$(basename "$f")" 'NR%4==2{n++; s+=length($0)} END{printf "%s\t%d\t%d\t%.1f\n", fn, n, s, (n? s/n:0)}' "$f"
    done
  } > reads/${SRR}_stats.txt
fi

fastqc -o reads reads/${SRR}*.fastq || true  # dont fail the script if FastQC prints warnings

# 8) Optional: clean up the tiny pilot subset
rm -rf pilot

echo "Done."
echo "Chosen run: $SRR (layout: $layout)"
echo "Pilot-measured bases/spot: ${bases_per_spot}"
echo "Target bases: ${target_bases} -> Spots downloaded: ${spots_needed}"
echo "Stats: reads/${SRR}_stats.txt ; FastQC HTML in reads/"

```
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

I searched for Ebola virus genomes sequenced by Nanopore and chose a random strain (SRR10769653).

```bash
# Set flags for ease of coding
set -xeuo pipefail

# ------------- User inputs (edit here) -------------
fastq=SRR10769653      # Nanopore run accession
genome_size=18959      # bp (Ebola ~18.959 kb)
coverage=10            # target depth (×10)
# ---------------------------------------------------

# ------ No edits beyond this point ------

# 1) Pilot: grab 100 reads to estimate real mean length
mkdir -p pilot reads
fastq-dump -X 100 --outdir pilot "$fastq"

# 2) Measure mean read length from pilot (single-end)
mean_len=$(awk 'NR%4==2{n++; s+=length($0)} END{printf "%.0f",(n?s/n:0)}' "pilot/${fastq}.fastq")
[ -z "$mean_len" ] || [ "$mean_len" -le 0 ] && { echo "Failed to estimate mean read length."; exit 1; }

# 3) Compute reads needed for ~10× (ceil division)
target_bases=$(( genome_size * coverage ))
num_reads=$(( (target_bases + mean_len - 1) / mean_len ))

echo "Mean read length ≈ ${mean_len} bp; target_bases=${target_bases}; downloading ~${num_reads} reads for ~${coverage}×"

# 4) Final download: only what you need (single FASTQ; no --split-files)
fastq-dump -X "$num_reads" --outdir reads "$fastq"

# 5) QC: stats + FastQC
if command -v seqkit >/dev/null 2>&1; then
  seqkit stats "reads/${fastq}.fastq" > "reads/${fastq}_stats.txt"
else
  # awk fallback
  awk 'NR%4==2{n++; s+=length($0)} END{
    printf "file\treads\ttotal_bases\tavg_len\n";
    printf "%s\t%d\t%d\t%.1f\n", FILENAME, n, s, (n? s/n:0)
  }' "reads/${fastq}.fastq" > "reads/${fastq}_stats.txt"
fi
fastqc -o reads "reads/${fastq}.fastq" || true

# 6) Clean up pilot
rm -rf pilot

# 7) Summary
echo "Done."
echo "Stats: reads/${fastq}_stats.txt"
echo "FastQC: reads/${fastq}_fastqc.html"

```
```bash
file                     format  type  num_seqs  sum_len  min_len  avg_len  max_len
reads/SRR10769653.fastq  FASTQ   DNA        300   77,650       35    258.8      695
```

<img width="1434" height="796" alt="image" src="https://github.com/user-attachments/assets/368236aa-f562-4720-81dd-61e4b967ff8d" />

Compared to the Illumina dataset, the Nanopore dataset shows much lower per-base quality.
