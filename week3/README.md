# Week 3: Visualising genomic data

## Downloading bacterial genome of interest 

Fasta file:
```bash
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/bacteria/current/fasta/bacteria_0_collection/klebsiella_pneumoniae_subsp_pneumoniae_mgh_78578_gca_000016305/dna/Klebsiella_pneumoniae_subsp_pneumoniae_mgh_78578_gca_000016305.ASM1630v1.dna_rm.toplevel.fa.gz
gunzip Klebsiella_pneumoniae_subsp_pneumoniae_mgh_78578_gca_000016305.ASM1630v1.dna_rm.toplevel.fa.gz
mv Klebsiella_pneumoniae_subsp_pneumoniae_mgh_78578_gca_000016305.ASM1630v1.dna_rm.toplevel.fa kleb.fa
```
gff file:
```bash
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/bacteria/current/gff3/bacteria_0_collection/klebsiella_pneumoniae_subsp_pneumoniae_mgh_78578_gca_000016305/Klebsiella_pneumoniae_subsp_pneumoniae_mgh_78578_gca_000016305.ASM1630v1.62.gff3.gz
gunzip Klebsiella_pneumoniae_subsp_pneumoniae_mgh_78578_gca_000016305.ASM1630v1.62.gff3.gz
mv Klebsiella_pneumoniae_subsp_pneumoniae_mgh_78578_gca_000016305.ASM1630v1.62.gff3 kleb.gff
```

_Note: Using shorter filenames (kleb.fa, kleb.gff) makes subsequent commands easier to manage and avoids long repetitive paths._

## How big is the genome?

```bash
 seqkit stats kleb.fa 
file     format  type  num_seqs    sum_len  min_len  avg_len    max_len
kleb.fa  FASTA   DNA          6  5,694,894    3,478  949,149  5,315,120
```

The total genome size is 5,694,894 bp, which includes the 5.3 Mb chromosome plus several smaller plasmids.

## How many features of each type does the GFF file contain?

```bash
grep -v "^#" kleb.gff | cut -f 3 | sort | uniq -c
```
As outlined by the output, the GFF file countains the following feature counts:

```bash
5185 CDS
  14 biological_region
   1 chromosome
5296 exon
5185 gene
5185 mRNA
 111 ncRNA_gene
  25 rRNA
   5 region
  86 tRNA
```

## From your GFF file, separate the intervals of type "gene" or "transcript" into a different file. Show the commands you used to do this.

```bash
grep -E '^#|\t(gene|mRNA|transcript)\t' kleb.gff > kleb_gene.gff
```

## Visualize the simplified GFF in IGV as a separate track. Compare the visualization of the original GFF with the simplified GFF.


<img width="1440" height="290" alt="1" src="https://github.com/user-attachments/assets/88a7f505-8937-4422-8cf6-004df84deb76" />

When loaded in IGV, the simplified file (kleb_gene.gff) shows a much cleaner view with only genes displayed, compared to the dense original file (kleb.gff) which includes CDS, exons, and other annotations. This makes it easier to see gene structure at a glance without the clutter of all sub-features.

## Zoom in to see the sequences, expand the view to show the translation table in IGV. Note how the translation table needs to be displayed in the correct orientation for it to make sense.

Forward strand (notice the arrow next to Sequence →):
<img width="1440" height="341" alt="2" src="https://github.com/user-attachments/assets/9c5b84b8-5999-47f7-9053-ae8267f1cf4d" />

Reverse strand (notice the arrow next to Sequence ←): 
<img width="1440" height="450" alt="3" src="https://github.com/user-attachments/assets/9a58a782-b48e-446c-87e3-46d60c97d1d7" />

If viewed in the wrong orientation, start and stop codons appear scrambled and are biologically meaningless.

## Visually verify that the first coding sequence of a gene starts with a start codon and that the last coding sequence of a gene ends with a stop codon.

I chose the gene purT. The first CDS begins with the start codon ATG (coding for methionine), shown in green as M in IGV.  

<img width="1440" height="322" alt="start" src="https://github.com/user-attachments/assets/b1a88822-6b14-41ce-a82a-692560bde4d8" />

The last CDS ends with a stop codon TAA, shown in red as **\***.
<img width="1440" height="352" alt="stop" src="https://github.com/user-attachments/assets/38e5867c-1db2-45c1-abef-ec034c6b2402" />
