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
