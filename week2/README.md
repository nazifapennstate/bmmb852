# Week 2: Demonstrate data analysis at UNIX command line

## Downloading the GFF3 file

You can download it directly from Ensembl using the following command:

```bash
wget https://ftp.ensembl.org/pub/current_gff3/apteryx_haastii/Apteryx_haastii.aptHaa1.115.gff3.gz
```

# Tell us a bit about the organism

The Great spotted kiwi (Apteryx haastii) is a large, nocturnal, flightless bird endemic to New Zealand.

# How many sequence regions (chromosomes) does the file contain?

```bash
zgrep "##sequence-region" Apteryx_haastii.aptHaa1.115.gff3.gz | wc -l
```
Output: 4183

# Does this match with the expectation for this organism?

Birds typically have ~30–40 chromosomes. The kiwi GFF3 lists 4183 sequence regions, which indicates that this assembly is still fragmented into thousands of scaffolds/contigs rather than a chromosome-level assembly.

# How many features does the file contain?

```bash
zgrep -v "^#" Apteryx_haastii.aptHaa1.115.gff3.gz | wc -l
```
Output: 870,888

# How many genes are listed for this organism?

```bash
zgrep -v "^#" Apteryx_haastii.aptHaa1.115.gff3.gz | cut -f3 | grep -w "gene" | wc -l
```

Output: 16,674

# Is there a feature type that you may have not heard about before? What is the feature and how is it defined?

Yes, I was not aware of the feature type biological_region, which is used in Ensembl to help annotate regulatory or functional DNA regions outside of protein-coding genes. Examples include CpG islands and predicted promoter regions.

# What are the top-ten most annotated feature types across the genome?

```bash
zgrep -v "^#" Apteryx_haastii.aptHaa1.115.gff3.gz | cut -f3 | sort | uniq -c | sort -nr | head -10
```

Output: 
```
337842 exon
329455 CDS
135966 biological_region
27823 mRNA
16674 gene
8251 five_prime_UTR
5247 three_prime_UTR
4183 region
2537 ncRNA_gene
2193 lnc_RNA
```

# Having analyzed this GFF file, does it seem like a complete and well-annotated organism?

The genome appears to be well-annotated but not fully assembled.

# Other insights 

The gene count (~16.6k), presence of multiple transcript types, and diversity of features (coding, non-coding, pseudogenes, biological region) indicate a thorough annotation. The very high number of sequence regions (4183) shows the genome is fragmented and not at chromosome resolution.

