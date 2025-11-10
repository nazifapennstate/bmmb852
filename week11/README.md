# Week 11: Establish the effects of variants

I ran my Makefile with `make batch`. That:

* pulled the Ebola reference (`GCF_000848505.1`) and indexed it,
* downloaded each of the samples' reads from SRA,
* aligned with `bwa mem`, sorted/indexed BAMs with `samtools`,
* called variants per sample with `bcftools`,
* and merged per-sample VCFs into `results/merged/merged.vcf.gz`.

I verified the merge with:

```bash
bcftools view -H results/merged/merged.vcf.gz | wc -l
bcftools stats results/merged/merged.vcf.gz | head
```

The merged file had **577 variant records**.

Then I annotated with snpEff:

```bash
make annotate_merged
```

This built a local snpEff DB from the FASTA+GFF and produced `results/merged/merged.snpeff.vcf.gz`.

To annotate **each** sample and get quick summaries:

```bash
awk -F, 'NR>1{printf "make effects SRR=%s SAMPLE=%s LAYOUT=%s\n",$1,$2,$3}' design.csv | bash
mkdir -p results/merged
cat results/*/*_effects_summary.txt > results/merged/all_effects_summary.txt
```

Then, I generated a table with **all** annotated variants across samples:

```bash
mkdir -p results/merged
out="results/merged/all_variants_annotated.tsv"
printf "Sample\tCHROM\tPOS\tREF\tALT\tQUAL\tEffect\tImpact\tGene\tAnnotation\n" > "$out"

while IFS=, read -r run sample layout; do
  [ "$run" = "Run" ] && continue
  vcf="results/$sample/${sample}.snpeff.vcf.gz"
  echo "Processing $sample..."
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%INFO/ANN\n' "$vcf" |
  awk -vS="$sample" -F'\t' '{
      split($6,a,"|"); effect=a[2]; impact=a[3]; gene=a[4];
      print S, $1, $2, $3, $4, $5, effect, impact, gene, $6
  }' OFS='\t' >> "$out"
done < design.csv
```

**Result:** `results/merged/all_variants_annotated.tsv`

---

## What the output shows:

<img width="1376" height="304" alt="image" src="https://github.com/user-attachments/assets/fa1f60e7-5479-419d-8ea4-3744d771e26e" />

## Column key (for `all_variants_annotated.tsv`)

| Column         | Description                                                                                                   | Example                                              | Interpretation                                                                                |
| :------------- | :------------------------------------------------------------------------------------------------------------ | :--------------------------------------------------- | :-------------------------------------------------------------------------------------------- |
| **Sample**     | Sample ID from `design.csv`. Identifies which isolate the variant came from.                                  | `EM096`                                              | Variant detected in sample EM096.                                                             |
| **CHROM**      | Reference contig or chromosome name.                                                                          | `NC_002549.1`                                        | All variants are aligned to the same Ebola reference genome.                                  |
| **POS**        | Genomic position (1-based).                                                                                   | `491`                                                | Variant occurs at nucleotide 491 of the reference genome.                                     |
| **REF**        | Reference base.                                                                                               | `A`                                                  | Base present in the reference genome.                                                         |
| **ALT**        | Alternate base in the sample.                                                                                 | `G`                                                  | Base found in the sample instead of reference.                                                |
| **QUAL**       | Phred-scaled quality score. Higher = greater confidence.                                                      | `225.4`                                              | High-confidence variant call.                                                                 |
| **Effect**     | Predicted functional effect from snpEff.                                                                      | `missense_variant`                                   | Describes the type of change (missense, synonymous, UTR, etc.).                               |
| **Impact**     | Severity category assigned by snpEff.                                                                         | `MODERATE`                                           | `MODIFIER` = non-coding, `LOW` = silent, `MODERATE` = amino-acid change, `HIGH` = disruptive. |
| **Gene**       | Gene affected by the variant.                                                                                 | `NP`                                                 | NP = nucleoprotein                            |
| **Annotation** | Full snpEff annotation (pipe-separated fields). Includes transcript, coding change, and protein substitution. | `A\|missense_variant\|MODERATE\|NP\|...p.Ile8Val...` | Provides the detailed predicted consequence on the gene/protein.                              |
