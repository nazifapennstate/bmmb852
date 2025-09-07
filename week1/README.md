# Week 1: Set up your system and demonstrate basic UNIX command line actions

# Check samtools version in the bioinfo environment

First, we need to activate our environment:

```
conda activate bioinfo
```


Now, check the samtools version:

```
samtools --version
```

Our output is:

```
samtools 1.22.1
```

# Show commands needed to create a nested directory structure

```
mkdir -p project/data/raw/fastq
```

# Show commands that create files in different directories

```
touch project/data/raw/sample1.txt
touch project/data/raw/fastq/sample_R1.fastq
```


# Show how to access these files using relative and absolute paths

Note that relative paths will depend on where your current pwd is:

```
cat project/data/raw/sample1.txt
cat project/data/raw/fastq/sample_R1.fastq
```

# Accessing using absolute paths

```
cat /Users/nazifa/Desktop/Fall-2025/Applied-Bioinformatics/bmmb852-week1/project/data/raw/sample1.txt
cat /Users/nazifa/Desktop/Fall-2025/Applied-Bioinformatics/bmmb852-week1/project/data/raw/fastq/sample_R1.fastq
```
