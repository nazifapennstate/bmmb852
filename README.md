# bmmb852-week1
Week 1 Assignment for BMMB852: Applied Bioinformatics 

Checking samtools version in the bioinfo environment

#Activate environment

conda activate bioinfo

#Check samtools version

samtools --version

#samtools 1.22.1

Answer to question: What version is your samtools command in the bioinfo environment? -samtools 1.22.1

Show commands needed to create a nested directory structure

mkdir -p project/data/raw/fastq

Show commands that create files in different directories

touch project/data/raw/sample1.txt
touch project/data/raw/fastq/sample_R1.fastq

Show how to access these files using relative and absolute paths

#Accessing using relative paths 

#This will depend on where your pwd is, mine is ~/Desktop/Fall-2025/Applied-Bioinformatics/bmmb852-week1 so I can do the following

cat project/data/raw/sample1.txt
cat project/data/raw/fastq/sample_R1.fastq

#Accessing using absolute paths

cat /Users/nazifachrf/Desktop/Fall-2025/Applied-Bioinformatics/bmmb852-week1/project/data/raw/sample1.txt
cat /Users/nazifachrf/Desktop/Fall-2025/Applied-Bioinformatics/bmmb852-week1/project/data/raw/fastq/sample_R1.fastq