# Week 11: Establish the effects of variants

I ran my Makefile pipeline with `make batch`. That:

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


## Column key (for `all_variants_annotated.tsv`)

* **Sample** — which sample the call came from (e.g., EM096).
* **CHROM** — contig name; here it’s the Ebola genome (`NC_002549.1`).
* **POS** — 1-based genomic position of the variant.
* **REF / ALT** — reference and alternate allele(s).
* **QUAL** — Phred-scaled variant quality from bcftools (higher → more confident call).
* **Effect** — snpEff’s top-ranked consequence for that allele (e.g., `5_prime_UTR_variant`, `synonymous_variant`, `missense_variant`).
* **Impact** — snpEff’s coarse severity bucket (`MODIFIER`, `LOW`, `MODERATE`, `HIGH`).
* **Gene** — the affected gene per annotation (e.g., `NP`, `VP35`, etc.).
* **Annotation** — the **full snpEff ANN field**: It includes the effect, impact, gene, transcript, coding change (e.g., `c.22A>G`) and protein change (e.g., `p.Ile8Val`).