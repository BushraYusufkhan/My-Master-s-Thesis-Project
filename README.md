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
mkdir fastqc_reports
fastqc *.fastq.gz -o ./fastqc_reports/
cd fastqc_reports/
mkdir multiqc_summary
multiqc ./ -o ./multiqc_summary/
```
Trimming of the reads:
```
mkdir -p cutadapt_trim_reads  # Create output folder if it doesn't exist

for R1 in *_R1_001.fastq.gz; do
    R2=${R1/_R1_001.fastq.gz/_R2_001.fastq.gz}                      
    OUT1=cutadapt_trim_reads/trimmed_${R1}
    OUT2=cutadapt_trim_reads/trimmed_${R2}                  

    echo "Processing $R1 and $R2 ..."

    cutadapt \
      -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
      -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
      --nextseq-trim=2 \
      -n 5 \
      -O 5 \
      -q 10,10 \
      -m 35 \
      -o $OUT1 \
      -p $OUT2 \
      $R1 $R2
done
```
For one of the samples, I used different value for -n and -q.
```
cutadapt \
  -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
  -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
  --nextseq-trim=2 \
  --trim-n \
  -n 7\
  -O 5 \
  -q 20,20 \
  -m 35 \
  -o trimmed_20-F314-Nla_S16_L001_R1_001.fastq.gz \
  -p trimmed_20-F314-Nla_S16_L001_R2_001.fastq.gz \
  20-F314-Nla_S16_L001_R1_001.fastq.gz 20-F314-Nla_S16_L001_R2_001.fastq.gz
```

