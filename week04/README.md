# Week 4: Obtain genomic data via accession numbers

I was assigned to Group 1, for which the paper was Gire et al. (2014), "Genomic surveillance elucidates Ebola virus origin and transmission during the 2014 outbreak," published in [_Science_](https://doi.org/10.1126/science.1259657).

## Downloading genome of interest 

Accession ID (GCF_000848505.1) was obtained by searching NCBI’s Genome/Assembly database for *Zaire ebolavirus* reference genomes. The record GCF_000848505.1 corresponds to the curated RefSeq assembly of the 1976 Yambuku-Mayinga isolate (NC_002549.1), which is widely used as the reference genome for Ebola virus.  

```bash
# Downloading fasta and gff file
datasets download genome accession GCF_000848505.1 \
  --include genome,gff3 \
  --filename GCF_000848505.1.zip

# Unzip into a temporary folder
unzip -o GCF_000848505.1.zip -d tmp_unzip

# Move & rename directly into the ebola folder
mv tmp_unzip/ncbi_dataset/data/*/*.fna  ebola.fasta
mv tmp_unzip/ncbi_dataset/data/*/genomic.gff ebola.gff

# Clean up the temporary folder and zip if you don’t need them
rm -rf tmp_unzip GCF_000848505.1.zip
```

## Genome visualisation using IGV 
<img width="1303" height="224" alt="Screenshot 2025-09-20 at 5 37 03 PM" src="https://github.com/user-attachments/assets/0f657f9d-0be8-4d5d-b4bf-f562e10e0a28" />

## How big is the genome?

```bash
grep -v "^>" ebola.fasta | tr -d '\n' | wc -c
```
Output:

```bash
18959
```

## How many features of each type does the GFF file contain?

```bash
awk '{print $3}' ebola.gff | sort | uniq -c
```
Output:
```bash
   6 
   1 1
  11 CDS
   1 annotwriter
   7 exon
   1 five_prime_UTR
   7 gene
   7 mRNA
   8 polyA_signal_sequence
   1 region
   8 regulatory_region
   4 sequence_feature
   1 three_prime_UTR
```
## What is the longest gene? What is its name and function?

```bash
awk '$3=="CDS"{print $5-$4+1, $0}' ebola.gff | sort -nr | head -1
```
Output:
```bash
6639 NC_002549.1	RefSeq	CDS	11581	18219	.	+	0	ID=cds-NP_066251.1;Parent=rna-ZEBOVgp7;Dbxref=GenBank:NP_066251.1,GeneID:911824;Name=NP_066251.1;gbkey=CDS;gene=L;locus_tag=ZEBOVgp7;product=RNA-dependent RNA polymerase;protein_id=NP_066251.1
(bioinfo) 
```
The longest gene is L, with a coding sequence of 6,639 nucleotides. It encodes the RNA-dependent RNA polymerase, the essential enzyme that copies the viral RNA genome and transcribes it into messenger RNAs.  

## Pick another gene and describe its name and function.

VP24 encodes a small viral protein of approximately 250 amino acids. It contributes to nucleocapsid stability and acts as an interferon antagonist by blocking STAT1 nuclear import, thereby suppressing the host antiviral response ([UniProt Q05322](https://www.uniprot.org/uniprotkb/Q05322/entry)).

## Looking at IGV, are the genomes closely packed or is there a lot of intergenomic space? Using IGV, estimate how much of the genome is covered by coding sequences.

From the IGV output, the genes appear tightly packed with only short intergenic regions. An estimated 90–95% of the genome is covered by coding sequences, consistent with the compact organization of negative-sense RNA viruses.


## Find alternative genome builds that could be used to perhaps answer a different question (find their accession numbers). Considering the focus of the paper, think about what other questions you could answer if you used a different genome build.

In addition to the RefSeq Mayinga 1976 reference (GCF_000848505.1) used here, another set of genomes is available in NCBI GenBank under accession codes ERX5245591–ERX5245598, MZ424849–MZ424862, and MZ605320–MZ605321. The analysis of these sequences was published in [Nature](https://doi.org/10.1038/s41586-021-03901-9). Using this build could enable new questions such as, how many mutations (SNPs) have accumulated since the 2014 outbreak and whether there are amino acid changes in known functional genes (such as VP24, and L) that might influence viral behavior.



