# My_Masters_Thesis_Project
### HPV integration revisited: Applying long-read sequencing data to study HPV integration in repeat-rich genomic regions
#### Download Human T2T Reference Genome

```bash
wget --content-disposition "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/GCF_009914755.1/download?include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GFF&include_annotation_type=RNA_FASTA&include_annotation_type=CDS_FASTA&include_annotation_type=PROT_FASTA&include_annotation_type=SEQUENCE_REPORT&hydrated=FULLY_HYDRATED"
```
```
unzip ncbi_dataset.zip
```

#### Download HPV Reference Genomes

HPV16
```bash
wget "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/GCF_000863945.3/download?include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GFF&include_annotation_type=RNA_FASTA&include_annotation_type=CDS_FASTA&include_annotation_type=PROT_FASTA&include_annotation_type=SEQUENCE_REPORT&hydrated=FULLY_HYDRATED" -O hpv16_dataset.zip
```
HPV18
```
curl -L -o genome_data.zip "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/GCF_000865665.1/download?include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GFF&include_annotation_type=RNA_FASTA&include_annotation_type=CDS_FASTA&include_annotation_type=PROT_FASTA&include_annotation_type=SEQUENCE_REPORT&hydrated=FULLY_HYDRATED"
```
### Datasets
#### Multi-omics mapping of human papillomavirus integration sites illuminates novel cervical cancer target genes (Iden at al)
Download the data from Iden et al (PRJNA640649)
```
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR120/055/SRR12056155/SRR12056155_subreads.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR120/056/SRR12056156/SRR12056156_subreads.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR120/057/SRR12056157/SRR12056157_subreads.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR120/058/SRR12056158/SRR12056158_subreads.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR120/059/SRR12056159/SRR12056159_subreads.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR120/060/SRR12056160/SRR12056160_subreads.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR120/061/SRR12056161/SRR12056161_subreads.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR120/062/SRR12056162/SRR12056162_subreads.fastq.gz

```
#### Proximity ligation‐based sequencing for the identification of human papillomavirus genomic integration sites in formalin‐fixed paraffin- embedded oropharyngeal squamous cell carcinomas (Demers et al., 2024) 
#### Quality control:

Before trimming:
```
module avail fastqc
module load bio/FastQC/0.12.1-Java-11
module avail multiqc
module load bio/MultiQC/1.22.3-foss-2023b
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar
```
Go to the working directory.
```
mkdir fastqc_reports
fastqc *.fastq.gz -o ./fastqc_reports/
cd fastqc_reports/
mkdir multiqc_summary
multiqc ./ -o ./multiqc_summary/
```
Trimming of the reads:
Go to the working directory.
```
mkdir trimmed_reads

for r1 in *_R1_001.fastq.gz; do
  r2="${r1/_R1_/_R2_}"
  sample="${r1%%_R1_001.fastq.gz}"

  echo "Processing $sample"

  java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads 8 -phred33 \
    "$r1" "$r2" \
    "trimmed_reads/${sample}_R1_paired_trimmed.fastq.gz" "trimmed_reads/${sample}_R1_unpaired.fastq.gz" \
    "trimmed_reads/${sample}_R2_paired_trimmed.fastq.gz" "trimmed_reads/${sample}_R2_unpaired.fastq.gz" \
    ILLUMINACLIP:$EBROOTTRIMMOMATIC/adapters/TruSeq3-PE.fa:2:30:10 \
    SLIDINGWINDOW:4:20 MINLEN:50
done

cd trimmed_reads/
mkdir unpaired_reads
mv *_unpaired.fastq.gz unpaired_reads/
mkdir after_cleaning_fastqc_reports
fastqc ./*_paired_trimmed.fastq.gz -o after_cleaning_fastqc_reports/
cd after_cleaning_fastqc_reports/
mkdir multiqc_trimmed_summary
multiqc ./ -o multiqc_trimmed_summary
```
Trim the partially trimmed adapters with Cutadapt.
Go to the Working directory.

Activate the Cutadapt environment.
conda activate cutadapt
```
mkdir -p ../cutadapt_cleaned_reads

for r1 in *_R1_paired_trimmed.fastq.gz; do
  r2="${r1/_R1_paired_trimmed.fastq.gz/_R2_paired_trimmed.fastq.gz}"
  sample=$(basename "$r1" _R1_paired_trimmed.fastq.gz)

  echo "Processing $sample with Cutadapt"

  cutadapt \
    -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
    -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    --minimum-length 50 --quality-cutoff 20 \
    -o ../cutadapt_cleaned_reads/${sample}_R1_final.fastq.gz \
    -p ../cutadapt_cleaned_reads/${sample}_R2_final.fastq.gz \
    "$r1" "$r2"
done
```

