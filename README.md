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
#### Mapping:
```
#!/bin/bash

REF="/scratch/bkhan1/hs1_HPV16/chm13v2_HPV16.fa"
THREADS=8
OUTDIR="/scratch/bkhan1/20240820_HPV_FFPE_TLC/combined_samples/Trimmed_reads_cutadapt/clean_reads_mapping"
READDIR="/scratch/bkhan1/20240820_HPV_FFPE_TLC/combined_samples/Trimmed_reads_cutadapt"

mkdir -p "$OUTDIR"

for R1 in "$READDIR"/*_R1_001.fastq.gz; do
    SAMPLE=$(basename "$R1" "_R1_001.fastq.gz")
    R2="$READDIR/${SAMPLE}_R2_001.fastq.gz"

    OUT_BAM="$OUTDIR/${SAMPLE}.sorted.bam"
    OUT_BAI="$OUTDIR/${SAMPLE}.sorted.bam.bai"

    if [[ -f "$OUT_BAM" && -f "$OUT_BAI" ]]; then
        echo "Sample $SAMPLE already processed, skipping."
        continue
    fi

    echo "Processing sample: $SAMPLE"

    bwa mem -M -t $THREADS -R "@RG\tID:$SAMPLE\tSM:$SAMPLE\tPL:ILLUMINA" "$REF" "$R1" "$R2" | \
    samtools view -b - | \
    samtools sort -@ $THREADS -o "$OUT_BAM"

    samtools index "$OUT_BAM"
done
```
#### Remove PCR duplicates:
```
# Path to Picard
PICARD="/projects/sw/eb/arch/zen4/software/picard/3.0.0-Java-17/picard.jar"

# Loop through all BAM files in the current folder
for bam in *.bam; do
    # Skip files that are already deduplicated
    if [[ $bam == *dedup.bam ]]; then
        continue
    fi

    # Define output names
    base=${bam%.bam}
    out="${base}.dedup.bam"
    metrics="${base}.metrics.txt"

    echo "Processing $bam ..."

    # Run Picard MarkDuplicates
    java -jar "$PICARD" MarkDuplicates \
        I="$bam" \
        O="$out" \
        M="$metrics" \
        REMOVE_DUPLICATES=true

    echo "Finished: $out"
    echo "Metrics saved to: $metrics"
    echo "----------------------------------------"
done

echo "All BAM files processed!"
```
#### Extract chimeric reads:

```
#!/bin/bash

# Create output folder if it doesn't exist
mkdir -p chimeric_reads

# Loop through all sorted BAM files in the current directory
for bam in *.sorted.dedup.bam; do
  # Extract base name (remove .sorted.dedup.bam)
  base=$(basename "$bam" .sorted.dedup.bam)

  echo "Processing $bam..."

  # Extract chimeric reads and write to BAM file
  samtools view -h "$bam" | \
  awk '
    /^@/ {print; next}
    /SA:Z:/ {
      main_chr = $3
      if (match($0, /SA:Z:([^,]+),/, m)) {
        sa_chr = m[1]
        if ((main_chr == "NC_001526.4" && sa_chr ~ /^chr/) ||
            (main_chr ~ /^chr/ && sa_chr == "NC_001526.4")) {
          print
        }
      }
    }
  ' | samtools view -b - > "chimeric_reads/${base}.chimeric.bam"

  # Index the new BAM file
  samtools index "chimeric_reads/${base}.chimeric.bam"
done
```
#### Find chimeric reads with Nla111 motif:
```
#!/bin/bash

# Script to detect NlaIII motifs (CATG / GTAC)
# Now also calculates distances between soft-clip boundary and motif site and combines unique read IDs

input_dir="."  # current folder
output_dir="10window"
mkdir -p "$output_dir"

for bam_file in "$input_dir"/*.bam; do
    name=$(basename "$bam_file" .bam)
    echo -e "\n============================"
    echo "Processing $name..."
    echo "============================"

    # === Step 1: Extract read types ===
    viral_human_sam="$output_dir/${name}_viralPrimary_humanSupp.sam"
    human_viral_sam="$output_dir/${name}_humanPrimary_viralSupp.sam"

    samtools view "$bam_file" | awk '$3 ~ /^NC_/ && $0 ~ /SA:Z:.*chr/' > "$viral_human_sam"
    samtools view "$bam_file" | awk '$3 ~ /^chr/ && $0 ~ /SA:Z:.*NC_/' > "$human_viral_sam"

    # === Step 2: Function to detect soft-clips + motifs ===
    process_reads () {
        local sam_file=$1
        local prefix=$2

        echo -e "\nExample CIGAR strings with soft-clipping from $prefix:"
        awk '{print $6}' "$sam_file" | grep -P '\d+S' | sort | uniq -c | sort -nr | head

        # --- Detect junction boundary from correct soft-clip side only ---
        awk '
        {
            read_id = $1;
            cigar = $6;
            seq = $10;

            # Find SA:Z tag
            sa_tag = "";
            for (i = 12; i <= NF; i++) {
                if ($i ~ /^SA:Z:/) {
                    sa_tag = substr($i, 6);
                    break;
                }
            }

            if (sa_tag == "") next;

            n = split(sa_tag, fields, ",");
            supp_pos = fields[2] + 0;
            supp_strand = fields[3];
            supp_cigar = fields[4];

            read_len = length(seq);

            # Determine which side aligns to supplementary segment
            if (match(cigar, /^[0-9]+S/)) {
                split(cigar, a, /[A-Z]/);
                left_clip = a[1];
            } else {
                left_clip = 0;
            }

            if (match(cigar, /[0-9]+S$/)) {
                match(cigar, /([0-9]+)S$/);
                right_clip = substr(cigar, RSTART, RLENGTH - 1);
            } else {
                right_clip = 0;
            }

            # Determine true junction side using strand logic
            if (supp_strand == "+") {
                # viral aligns to right of read
                if (right_clip > 0) {
                    alignment_end = read_len - right_clip;
                    print read_id "\t" alignment_end;
                }
            } else {
                # viral aligns to left of read
                if (left_clip > 0) {
                    print read_id "\t" left_clip;
                }
            }
        }
        ' "$sam_file" > "${sam_file%.sam}_softclip_boundaries.txt"

        # --- Motif positions (CATG/GTAC in read) ---
        cut -f1,10 "$sam_file" | perl -ne '
            chomp;
            my ($id, $seq) = split(/\t/);
            while ($seq =~ /CATG/g) {
                print "$id\t" . (pos($seq) - 3) . "\n";
            }
            while ($seq =~ /GTAC/g) {
                print "$id\t" . (pos($seq) - 3) . "\n";
            }
        ' > "${sam_file%.sam}_motif_positions.txt"

        # --- Match soft-clip boundary with motif positions (within 10bp) ---
        join -t $'\t' -1 1 -2 1 \
            <(sort -k1,1 "${sam_file%.sam}_softclip_boundaries.txt") \
            <(sort -k1,1 "${sam_file%.sam}_motif_positions.txt") \
        | awk -F'\t' '
            function abs(v) { return v < 0 ? -v : v }
            {
                clip = $2;
                motif = $3;
                if (abs(clip - motif) <= 10) {
                    print $0;
                }
            }
        ' > "${sam_file%.sam}_IntSites_NlaIII_sites.txt"

        echo -e "\n[$prefix] Reads matching motif within 5bp of soft-clip:"
        wc -l "${sam_file%.sam}_IntSites_NlaIII_sites.txt"

        # --- Count unique read IDs ---
        cut -f1 "${sam_file%.sam}_IntSites_NlaIII_sites.txt" | sort | uniq > "${sam_file%.sam}_IntSites_NlaIII_uniqueReads.txt"
        echo -e "[$prefix] Total unique reads at NlaIII junctions:"
        wc -l "${sam_file%.sam}_IntSites_NlaIII_uniqueReads.txt"

        # --- Also save softclip-motif distances for plotting ---
        join -t $'\t' -1 1 -2 1 \
            <(sort -k1,1 "${sam_file%.sam}_softclip_boundaries.txt") \
            <(sort -k1,1 "${sam_file%.sam}_motif_positions.txt") \
        | awk -F'\t' '
            {
                clip = $2;
                motif = $3;
                dist = clip - motif;
                print $1 "\t" clip "\t" motif "\t" dist;
            }
        ' > "${sam_file%.sam}_IntSites_NlaIII_distances.txt"

        echo -e "[$prefix] Saved distances to: ${sam_file%.sam}_IntSites_NlaIII_distances.txt"
    }

    # Process both chimeric directions
    process_reads "$viral_human_sam" "${name}_viralPrimary_humanSupp"
    process_reads "$human_viral_sam" "${name}_humanPrimary_viralSupp"

    # === Combine unique reads from both directions ===
    echo -e "\n[Combining unique reads from both directions...]"
    cat "${viral_human_sam%.sam}_IntSites_NlaIII_uniqueReads.txt" "${human_viral_sam%.sam}_IntSites_NlaIII_uniqueReads.txt" \
        | sort | uniq > "$output_dir/${name}_combined_uniqueReads.txt"

    echo -e "[Combined] Total unique reads at junctions (viral→human + human→viral):"
    wc -l "$output_dir/${name}_combined_uniqueReads.txt"

done
```
#### Exclude the reads which have Nla111 restriction motif within a window of 10bp upstream and downstream of the split junction:
```
#!/bin/bash

READ_LIST_DIR="/scratch/bkhan1/20240820_HPV_FFPE_TLC/combined_folder/samples-cut-with-nla/cutadapt_trim_reads/clean_reads_mapping/ded
up-files/chimeric_reads/10window/reads-to-exclude"
TRUE_INTEGRATION_DIR="Filtered-bam-files"
NLA_MOTIF_DIR="Bam-files-of-reads-which-have-nla-motif-within-window10"
mkdir -p "$TRUE_INTEGRATION_DIR" "$NLA_MOTIF_DIR"

for BAM_FILE in *.bam; do
    BASENAME="${BAM_FILE%.bam}"

    # Find read list by partial match
    READ_LIST=$(ls "$READ_LIST_DIR"/*"$BASENAME"* 2>/dev/null | head -n 1)

    if [[ -z "$READ_LIST" ]]; then
        echo "No matching read list found for $BASENAME — skipping."
        continue
    fi

    echo "Processing $BASENAME with read list $(basename "$READ_LIST")"

    # 1. BAM without NlaIII motif reads
    samtools view -h "$BAM_FILE" \
    | grep -vFf "$READ_LIST" \
    | samtools view -b -o "${TRUE_INTEGRATION_DIR}/${BASENAME}.true_integration.bam"

    # 2. BAM with NlaIII motif reads
    samtools view -h "$BAM_FILE" | \
    awk -v read_list="$READ_LIST" 'BEGIN {
        while ((getline line < read_list) > 0)
            reads[line] = 1
    }
    (/^@/) { print; next }
    { if ($1 in reads) print }
    ' | samtools view -bS - > "${NLA_MOTIF_DIR}/${BASENAME}.integration_with_nla_motif.bam"

done
```
#### Sort the reads into clusters to determine human HPV breakpoints:
```
#!/bin/bash

# Directory containing BAM files
BAM_DIR="./"   # change if your BAM files are elsewhere
OUT_DIR="./processed_samples"

mkdir -p "$OUT_DIR"

# Summary file for all samples
SUMMARY="$OUT_DIR/all_samples_cluster_counts.tsv"
echo -e "Sample\tClusterID\tReadCount" > "$SUMMARY"

# Loop over BAM files
for bam in "$BAM_DIR"/*.bam; do
    sample=$(basename "$bam" .bam)
    sample_dir="$OUT_DIR/$sample"
    mkdir -p "$sample_dir"

    echo "Processing sample: $sample"

    # Extract human reads (everything NOT NC_001526.4)
    samtools view -h "$bam" | awk '$3 != "NC_001526.4" || $1 ~ /^@/' | \
        samtools view -b -o "$sample_dir/${sample}_human.bam"

    # Extract viral reads (only NC_001526.4)
    samtools view -b "$bam" "NC_001526.4" > "$sample_dir/${sample}_viral.bam"

    # Convert BAM to BED
    bedtools bamtobed -i "$sample_dir/${sample}_human.bam" > "$sample_dir/${sample}_human.bed"
    bedtools bamtobed -i "$sample_dir/${sample}_viral.bam" > "$sample_dir/${sample}_viral.bed"

    # Sort BED by read name (4th column)
    sort -k4,4 "$sample_dir/${sample}_human.bed" > "$sample_dir/${sample}_human_sorted.bed"
    sort -k4,4 "$sample_dir/${sample}_viral.bed" > "$sample_dir/${sample}_viral_sorted.bed"

    # Optional join human and viral reads by read name
    join -1 4 -2 4 <(sort -k4,4 "$sample_dir/${sample}_human_sorted.bed") \
                    <(sort -k4,4 "$sample_dir/${sample}_viral_sorted.bed") \
                    > "$sample_dir/${sample}_human_viral_joined.bed"

    # Cluster human reads
    bedtools cluster -i "$sample_dir/${sample}_human.bed" > "$sample_dir/${sample}_human_clustered.tsv"

    # Count reads per cluster and append to summary
    awk -F'\t' -v sample="$sample" '{counts[$NF]++} END {for (c in counts) print sample"\t"c"\t"counts[c]}' \
        "$sample_dir/${sample}_human_clustered.tsv" >> "$SUMMARY"

    # Split clusters into separate files
    awk -F'\t' -v outdir="$sample_dir" -v sample="$sample" '{ out = outdir"/"sample"_cluster_" $NF ".tsv"; print > out }' \
        "$sample_dir/${sample}_human_clustered.tsv"

    echo "Finished processing $sample"
done

echo "All BAM files processed. Summary saved to $SUMMARY"
```
I was facing an issue with short reads; sometimes the reads had identical start and end positions. I did not know what those were. I had already removed PCR duplicates. In the original study, they said something about low-complexity hairpin regions. And I assumed the reads which had the same start and end would be these low complexity hairpin regions. And to get rid of those, I counted them as one; those reads were distorting the true signals.
```
base_dir="/scratch/bkhan1/20240820_HPV_FFPE_TLC/combined_folder/samples-cut-with-nla/cutadapt_trim_reads/clean_reads_mapping/dedup-fi
les/chimeric_reads/filtered-bam-files/processed_samples"

cd "$base_dir" || exit

# Output summary files
cluster_summary="3dedup_cluster_summary.tsv"
total_summary="3dedup_total_read_counts.tsv"

# Headers
echo -e "Sample\tCluster_File\tBefore_Dedup_Reads\tAfter_Dedup_Reads\tCount_Entries" > "$cluster_summary"
echo -e "Sample\tTotal_Reads" > "$total_summary"

# Loop through each sample folder
for sample_dir in trimmed_*; do
    if [ -d "$sample_dir" ]; then
        echo "------------------------------------------"
        echo " Processing sample: $sample_dir"
        echo "------------------------------------------"

        cd "$sample_dir" || continue
        mkdir -p dedup counts

        # Process all cluster files
        for file in *_cluster_*.tsv; do
            [ -s "$file" ] || { echo "Skipping empty file: $file"; continue; }

            dedup_file="dedup/${file%.tsv}_dedup.tsv"
            counts_file="dedup/${file%.tsv}_counts.tsv"

            # Deduplicate by first 3 columns
            sort -t$'\t' -k1,1 -k2,2n -k3,3n "$file" | awk '!seen[$1"\t"$2"\t"$3]++' > "$dedup_file"

            # Count reads per cluster
            awk -F'\t' '{key=$1"\t"$2"\t"$3; counts[key]++} END {
                for(k in counts) print k"\t"counts[k]
            }' "$dedup_file" > "$counts_file"

            # Before/After dedup read counts
            before_dedup=$(wc -l < "$file")
            after_dedup=$(wc -l < "$dedup_file")
            num_counts=$(wc -l < "$counts_file")

            # Append summary info
            echo -e "${sample_dir}\t${file}\t${before_dedup}\t${after_dedup}\t${num_counts}" >> "$base_dir/$cluster_summary"
        done

        # Safe version (avoids 'Argument list too long' error)
        total_reads=$(find . -maxdepth 1 -name "*_cluster_*.tsv" -type f -exec cat {} + | awk 'NF >= 3 {reads++} END{print reads}')
        echo -e "${sample_dir}\t${total_reads}" >> "$base_dir/$total_summary"

        cd "$base_dir" || exit
    fi
done

echo "Processing complete!"
echo "Summaries saved to:"
echo "   - $cluster_summary"
echo "   - $total_summary"
```
```
#!/bin/bash

# Input file
DEDUP_SUMMARY="3dedup_cluster_summary.tsv"

# Output files
OUT_1P="3dedup_clusters_1percent_min3reads.tsv"
OUT_01P="3dedup_clusters_0.01percent_min3reads.tsv"

# Remove old outputs
rm -f "$OUT_1P" "$OUT_01P"

# Headers (added Before_Dedup_Reads)
echo -e "Sample\tCluster_ID\tBefore_Dedup_Reads\tAfter_Dedup_Reads\tPercentage" > "$OUT_1P"
echo -e "Sample\tCluster_ID\tBefore_Dedup_Reads\tAfter_Dedup_Reads\tPercentage" > "$OUT_01P"

# Iterate over each unique sample
samples=$(awk 'NR>1 {print $1}' "$DEDUP_SUMMARY" | sort -u)

for sample in $samples; do
    # Clean up sample name
    clean_sample=$(echo "$sample" | sed -E 's/^trimmed_//; s/(_L001)?\.sorted\.dedup\.chimeric(\.true_integration)?//')

    # Total reads (sum of After_Dedup_Reads = column 4)
    total_reads=$(awk -v s="$sample" -F'\t' '$1==s {sum+=$4} END{print sum}' "$DEDUP_SUMMARY")

    # Skip if no reads
    [ -z "$total_reads" ] && continue

    # Extract and calculate per cluster
    awk -v s="$sample" -v total="$total_reads" -F'\t' '
        $1==s && $4>2 {
            cluster=$2
            sub(/.*_cluster_/, "", cluster)
            sub(/\.tsv$/, "", cluster)
            before=$3
            after=$4
            perc=(after/total)*100
            print s "\t" cluster "\t" before "\t" after "\t" perc
        }' "$DEDUP_SUMMARY" |
    while read -r s cluster before after perc; do
        # Write to 1% threshold file
        if (( $(echo "$perc >= 1" | bc -l) )); then
            echo -e "$clean_sample\t$cluster\t$before\t$after\t$perc" >> "$OUT_1P"
        fi
        # Write to 0.1% threshold file
        if (( $(echo "$perc >= 0.01" | bc -l) )); then
            echo -e "$clean_sample\t$cluster\t$before\t$after\t$perc" >> "$OUT_01P"
        fi
    done
done

echo "Done. Generated:"
echo "   - $OUT_1P"
echo "   - $OUT_01P"
```
