# ============================================================
# Makefile for KRAS region comparison across sequencing platforms
# ============================================================

# KRAS region from assignment
REGION = chr12:25205246-25250936

# Reference genome (GRCh38 GIABv3)
REF_URL = https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta.gz
REF = refs/GRCh38.fasta

# -----------------------
# Platform BAM URLs
# -----------------------

# Illumina N-D sample (118x, GRCh38-GIABv3)
ILLUMINA_BAM_URL = https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/NYGC_Illumina-WGS_20231023/HG008-N-D_Illumina_118x_GRCh38-GIABv3.bam

# Element AVITI N-D sample (77x, GRCh38-GIABv3)
ELEMENT_BAM_URL = https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/Element-AVITI-20241216/HG008-N-D_Element-StdInsert_77x_GRCh38-GIABv3.bam

# PacBio Revio N-P sample (35x, GRCh38-GIABv3)
LONGREAD_BAM_URL = https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/PacBio_Revio_20240125/HG008-N-P_PacBio-HiFi-Revio_20240125_35x_GRCh38-GIABv3.bam

# -----------------------
# Directories
# -----------------------
all: illumina element longread

refs:
	mkdir -p refs

bam:
	mkdir -p bam

# -----------------------
# Reference genome
# -----------------------
$(REF): | refs
	curl -L $(REF_URL) > $(REF).gz
	gunzip -f $(REF).gz
	samtools faidx $(REF)

# -----------------------
# Illumina
# -----------------------
illumina: $(REF) | bam
	samtools view -b $(ILLUMINA_BAM_URL) $(REGION) > bam/kras_illumina.bam
	samtools index bam/kras_illumina.bam
	samtools flagstat bam/kras_illumina.bam > bam/kras_illumina.flagstat.txt
	samtools depth -r $(REGION) bam/kras_illumina.bam > bam/kras_illumina.depth.txt

# -----------------------
# Element AVITI
# -----------------------
element: $(REF) | bam
	samtools view -b $(ELEMENT_BAM_URL) $(REGION) > bam/kras_element.bam
	samtools index bam/kras_element.bam
	samtools flagstat bam/kras_element.bam > bam/kras_element.flagstat.txt
	samtools depth -r $(REGION) bam/kras_element.bam > bam/kras_element.depth.txt

# -----------------------
# PacBio Revio (HiFi)
# -----------------------
longread: $(REF) | bam
	samtools view -b $(LONGREAD_BAM_URL) $(REGION) > bam/kras_longread.bam
	samtools index bam/kras_longread.bam
	samtools flagstat bam/kras_longread.bam > bam/kras_longread.flagstat.txt
	samtools depth -r $(REGION) bam/kras_longread.bam > bam/kras_longread.depth.txt

# -----------------------
# Cleanup
# -----------------------
clean:
	rm -rf refs bam