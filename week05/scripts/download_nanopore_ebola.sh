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
