# Week 7: Write a reusable alignment Makefile

This week’s goal was to generalize the Week 6 Makefile so that it can process both Illumina and Oxford Nanopore (ONT) sequencing runs for the Ebola virus reference genome. The pipeline automatically downloads data, aligns reads using the appropriate mapper, and generates coverage tracks for visualization.

| Target          | Purpose                                                                                                    |
| --------------- | ---------------------------------------------------------------------------------------------------------- |
| `all`           | Runs the full pipeline from genome download through BigWig generation                                      |
| `$(GENOME)`     | Downloads and extracts the Ebola reference genome (RefSeq GCF_000848505.1)                                 |
| `fastq`         | Downloads SRA reads (`fastq-dump` or `fasterq-dump`) — supports paired-end (Illumina) or single-end (ONT)  |
| `qc`            | Runs FastQC on reads (if available on PATH)                                                                |
| `index` / `mmi` | Builds genome indexes for BWA (Illumina) and minimap2 (ONT)                                                |
| `align`         | Chooses the correct aligner: `bwa mem` for Illumina or `minimap2 -ax map-ont` for ONT → sorted BAM + index |
| `stats`         | Calculates read-mapping statistics (`samtools flagstat`) and mean coverage (`samtools depth`)              |
| `bigwig`        | Converts BAM coverage to bedGraph and BigWig for IGV visualization                                         |
| `clean`         | Removes all generated files and directories                                                                |

After activating the `bioinfo` Conda environment, I ran pipeline with the following SRR:

```bash
# Illumina paired-end
make all SRR=SRR1734993 LAYOUT=PE PLATFORM=illumina SUBSET=

# Nanopore single-end
make all SRR=SRR8959866 LAYOUT=SE PLATFORM=ont SUBSET=
```

_Note: In previous assignments, I used `SRR10769653`, which was a **metagenomic GridION ONT run** with spiked primer enrichment.
However, it showed extremely poor coverage in IGV, so I replaced it with **SRR8959866**, another ONT run with better read depth.
For this week, I ran the pipeline **without specifying a subset** to see coverage beyond 10×, though the Makefile supports downsampling via `SUBSET=...` if needed (e.g., `SUBSET=100000` for ~10× coverage)._

All results are stored in structured subdirectories:

```
genome/           → reference FASTA + index files  
reads/            → raw FASTQ files + FastQC reports  
results/<SRR>/    → BAM + BAI + stats + BigWig  
```

### Changes to the Makefile this week - adding BigWig and using two aligners

BAM files show individual alignments, but BigWig files provide a **continuous coverage profile** across the genome.
Creating them with `bedtools genomecov` + `bedGraphToBigWig` allows quick visual comparisons of coverage depth between datasets in IGV.

Illumina short reads are best handled by **BWA-MEM**, which efficiently maps paired, high-accuracy reads.
Oxford Nanopore reads are long, error-prone, and single-end, requiring **minimap2** with the `map-ont` preset optimized for long-read alignments.

## Briefly describe the differences between the alignment in both files.

## Briefly compare the statistics for the two BAM files.

I compared the resulting alignment statistics with this:

```
for S in SRR1734993 SRR8959866; do
  ST=results/$S/${S}_alignment_stats.txt
  total=$(awk '/in total/ {print $1}' "$ST")
  primary=$(awk '/primary$/ {print $1}' "$ST")
  mapped=$(awk '/ mapped \(/ && !/supplementary/ {gsub(/[()%]/,"",$5); print $5}' "$ST")
  proper=$(awk '/properly paired/ {print $1}' "$ST")
  mean=$(awk '/Average coverage/ {print $4}' "$ST")
  echo -e "$S\t$total\t$primary\t$mapped\t$proper\t$mean" >> results/comparison_summary.tsv
done

column -t results/comparison_summary.tsv
```

| SRR ID         | Platform | Layout     | Total Reads | Primary Alignments | Mapped % | Properly Paired | Mean Coverage (×) |
| :------------- | :------- | :--------- | ----------: | -----------------: | -------: | --------------: | ----------------: |
| **SRR1734993** | Illumina | Paired-end |   1,090,855 |          1,090,342 |  ~1.15 % |             Yes |            56.3 × |
| **SRR8959866** | ONT      | Single-end |   1,205,102 |          1,122,688 |  ~86.6 % |             N/A |        26,422.5 × |


## How many primary alignments does each of your BAM files contain?

SRR1734993 (Illumina) has 1,090,342 and SRR8959866 (ONT) has 1,122,688.

## What coordinate has the largest observed coverage (hint samtools depth)

I used samtools depth to identify the position with the highest coverage in each BAM file:

```
for S in SRR1734993 SRR8959866; do
  BAM=results/$S/${S}_sorted.bam
  samtools depth "$BAM" | sort -k3,3nr | head -1 | awk -v s=$S '{printf "%s\t%s\t%s\t%s\n", s, $1, $2, $3}'
done
```
Results:
```
SRR1734993      NC_002549.1     18425   316
SRR8959866      NC_002549.1     8949    130087
```
The Illumina run (SRR1734993) shows its maximum coverage (~316×) at 18425 bp.

The Nanopore run (SRR8959866) shows its maximum coverage (~130k×) at position 8,949 bp.

## Select a gene of interest. How many alignments on a forward strand cover the gene?

Initially, I triedto identify and quantify gene coverage visually in IGV, but the alignments were too dense and strand separation was not easily interpretable. Instead, I examined the annotation file (genome/ebola.gff) using the less command and selected the first coding gene, NP (nucleoprotein; coordinates 56–3026 on NC_002549.1, forward strand).

To quantify forward-strand coverage over this gene, I ran the following command:

```
for S in SRR1734993 SRR8959866; do
  BAM=results/$S/${S}_sorted.bam
  REGION="NC_002549.1:56-3026"
  forward=$(samtools view -c -F 16 "$BAM" "$REGION")
  echo -e "$S\tForward_reads_over_NP\t$forward"
done
```

This produced:
```
SRR ID	Platform	Forward Reads Covering NP
SRR1734993	Illumina (PE)	897
SRR8959866	ONT (SE)	117,142
```

The ONT dataset showed far greater read depth over NP, consistent with the much higher overall coverage and single-end strand bias observed in Nanopore amplicon-based sequencing.