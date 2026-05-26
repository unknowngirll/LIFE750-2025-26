#!/bin/bash
# LIFE750 Cycle 2 
# All the commands I used for this analysis
# Student ID: 201930340


# Part 1: Variant calling on the mutant Gene X allele

# Activate the conda environment with all the tools
conda activate life750

# Move to the project directory
cd ~/LIFE750_project

# Make folders for the outputs
mkdir variant_calling
mkdir quality


#Quality control 
# Run FastQC on the four fastq files
fastqc --outdir quality --noextract --nogroup --format fastq *.fastq


# Index the Gene X reference 

# Index the reference for bwa
bwa index GeneX_reference.fa

# Index the reference for samtools (needed for IGV later)
samtools faidx GeneX_reference.fa


#  Map the normal sample to the Gene X reference

# Map with bwa mem, convert to BAM, sort, all piped together
bwa mem -t 2 -R '@RG\tID:normal\tSM:normal\tPL:ILLUMINA\tLB:WGS' \
  GeneX_reference.fa normal_GeneX_R1.fastq normal_GeneX_R2.fastq \
  | samtools view -b -@ 2 \
  | samtools sort -@ 2 -o variant_calling/normal_sorted.bam -

# Index the sorted BAM
samtools index variant_calling/normal_sorted.bam


# for mutant(same steps)

bwa mem -t 2 -R '@RG\tID:mutant\tSM:mutant\tPL:ILLUMINA\tLB:WGS' \
  GeneX_reference.fa mutant_GeneX_R1.fastq mutant_GeneX_R2.fastq \
  | samtools view -b -@ 2 \
  | samtools sort -@ 2 -o variant_calling/mutant_sorted.bam -

samtools index variant_calling/mutant_sorted.bam


# Check that mapping worked 

# Summary stats for both BAM files
samtools flagstat variant_calling/normal_sorted.bam
samtools flagstat variant_calling/mutant_sorted.bam


#Variant calling

# Call variants in the mutant sample with FreeBayes
# Using --ploidy 1 because the plasmid is haploid
freebayes --ploidy 1 --min-mapping-quality 30 --min-base-quality 25 \
  --fasta-reference GeneX_reference.fa \
  variant_calling/mutant_sorted.bam > variant_calling/mutant_variants.vcf

# Call variants in the normal sample (as a control - should be zero)
freebayes --ploidy 1 --min-mapping-quality 30 --min-base-quality 25 \
  --fasta-reference GeneX_reference.fa \
  variant_calling/normal_sorted.bam > variant_calling/normal_variants.vcf

# Count variants in each
grep -v "^#" variant_calling/normal_variants.vcf | wc -l
grep -v "^#" variant_calling/mutant_variants.vcf | wc -l


#Filter low-confidence variants
# Keep only variants with good quality, depth and allele frequency
bcftools view -e 'QUAL < 10 || INFO/DP < 10 || AF = 0' \
  -O v -o variant_calling/mutant_variants_filtered.vcf \
  variant_calling/mutant_variants.vcf

# How many passed the filter?
grep -v "^#" variant_calling/mutant_variants_filtered.vcf | wc -l


#  Build a consensus mutant Gene X sequence-

# bcftools consensus needs a bgzipped, indexed VCF
bgzip variant_calling/mutant_variants_filtered.vcf
bcftools index variant_calling/mutant_variants_filtered.vcf.gz

# Generate the consensus FASTA
bcftools consensus -f GeneX_reference.fa \
  -o variant_calling/mutant_GeneX_consensus.fa \
  variant_calling/mutant_variants_filtered.vcf.gz

# Rename the header so it doesn't say GeneX_reference
sed -i 's/>GeneX_reference/>GeneX_mutant/' variant_calling/mutant_GeneX_consensus.fa


#Inspect the sequence around each variant

# Show the reference and mutant base around position 413
samtools faidx GeneX_reference.fa GeneX_reference:408-418
samtools faidx variant_calling/mutant_GeneX_consensus.fa GeneX_mutant:408-418

# Same for position 573
samtools faidx GeneX_reference.fa GeneX_reference:568-578
samtools faidx variant_calling/mutant_GeneX_consensus.fa GeneX_mutant:568-578

# Same for position 617
samtools faidx GeneX_reference.fa GeneX_reference:612-622
samtools faidx variant_calling/mutant_GeneX_consensus.fa GeneX_mutant:612-622


#Annotate the variants with bcftools csq

# This adds the BCSQ tag (synonymous / missense / stop_gained)
bcftools csq --fasta-ref GeneX_reference.fa \
  --gff-annot GeneX.gff3 \
  variant_calling/mutant_variants_filtered.vcf.gz \
  > variant_calling/mutant_variants_annotated.vcf

# View the annotated variants
grep -v "^#" variant_calling/mutant_variants_annotated.vcf

# Summarise the consequences
egrep --only-matching "BCSQ=.+?\|" variant_calling/mutant_variants_annotated.vcf \
  | sort | uniq -c




#Part 2
#checking the data before Rstudio
cat sample_metadata.tsv
head -5 gene_counts.tsv
wc -l gene_counts.tsv
head -1 gene_counts.tsv | awk '{print NF}'


#part 3

cd ~/LIFE750_project
mkdir -p cutrun_analysis

#Prepare gene annotation file
# Extract only gene features from the main GFF annotation
awk '$3 == "gene"' genes.gff3 > genes_only.gff3

# Check how many genes are in the annotation
wc -l genes_only.gff3
# Inspect the first few lines to confirm the format
head genes_only.gff3

#Find genes bound by Gene X using bedtools intersect
# Use bedtools intersect to find overlap between gene annotations and Gene X binding sites
# -a is for  gene annotations
# -b is for binding sites
bedtools intersect -a genes_only.gff3 -b gene_x_binding_sites.bed -wa > genes_bound_by_geneX.gff3

# Count how many gene-binding site overlaps were found
wc -l genes_bound_by_geneX.gff3


# Inspect the first few overlapping genes
head genes_bound_by_geneX.gff3


# Extract unique Gene IDs from the GFF attribute column
# grep with Perl regex to capture text after "ID=" until ";"
grep -oP 'ID=\K[^;]+' genes_bound_by_geneX.gff3 | sort -u > geneX_target_genes.txt


# Count unique target genes
wc -l geneX_target_genes.txt

# View the target genes
head geneX_target_genes.txt

#Compare Gene X targets with differentially expressed genes

# Extract Gene IDs from significant_genes.tsv (skip the header line)
tail -n +2 cutrun_analysis/significant_genes.tsv | cut -f1 > cutrun_analysis/de_genes.txt

# Check the DE genes file
wc -l cutrun_analysis/de_genes.txt
head cutrun_analysis/de_genes.txt

# in RStudio on Windows the created files have some extra endings which I got error running 'comm'
# Convert them to Unix LF line endings.
sed -i 's/\r$//' cutrun_analysis/upregulated_genes.txt
sed -i 's/\r$//' cutrun_analysis/downregulated_genes.txt
sed -i 's/\r$//' cutrun_analysis/de_genes.txt

head -3 cutrun_analysis/upregulated_genes.txt | cat -A

# Find overlap between Gene X targets and DE genes

# 'comm -12' shows lines common to both sorted files
comm -12 <(sort geneX_target_genes.txt) <(sort cutrun_analysis/de_genes.txt) > cutrun_analysis/overlap_target_DE.txt

# View the overlap genes
wc -l cutrun_analysis/overlap_target_DE.txt
cat cutrun_analysis/overlap_target_DE.txt


# Split overlap into UP- and DOWN-regulated genes
# View the genes that are both Gene X targets AND differentially expressed
cat cutrun_analysis/overlap_target_DE.txt

# Find targets that are UP-regulated in mutant

comm -12 <(sort geneX_target_genes.txt) <(sort cutrun_analysis/upregulated_genes.txt) > cutrun_analysis/overlap_target_UP.txt
wc -l cutrun_analysis/overlap_target_UP.txt
cat cutrun_analysis/overlap_target_UP.txt

# Find targets that are DOWN-regulated in mutant

comm -12 <(sort geneX_target_genes.txt) <(sort cutrun_analysis/downregulated_genes.txt) > cutrun_analysis/overlap_target_DOWN.txt
wc -l cutrun_analysis/overlap_target_DOWN.txt
cat cutrun_analysis/overlap_target_DOWN.txt



# Numbers needed for enrichment calculation:
# Expected overlap = (targets x DE genes) / total genes
wc -l genes_only.gff3
wc -l geneX_target_genes.txt
wc -l cutrun_analysis/de_genes.txt
wc -l cutrun_analysis/overlap_target_DE.txt

#Compare real peaks vs decoy peaks

# Separate real peaks from decoy peaks
grep "^chr1" gene_x_binding_sites.bed | grep -v "decoy" > real_peaks.bed
grep "decoy" gene_x_binding_sites.bed > decoy_peaks.bed

wc -l real_peaks.bed
wc -l decoy_peaks.bed


# Find genes near real peaks only
bedtools intersect -a genes_only.gff3 -b real_peaks.bed -wa | grep -oP 'ID=\K[^;]+' | sort -u > geneX_real_targets.txt
wc -l geneX_real_targets.txt

# Find genes near decoy peaks only
bedtools intersect -a genes_only.gff3 -b decoy_peaks.bed -wa | grep -oP 'ID=\K[^;]+' | sort -u > geneX_decoy_targets.txt
wc -l geneX_decoy_targets.txt

# Fix line endings for these new files
sed -i 's/\r$//' geneX_real_targets.txt
sed -i 's/\r$//' geneX_decoy_targets.txt



# Count how many REAL targets are UP-regulated

# Real targets that are UP-regulated 
comm -12 <(sort geneX_real_targets.txt) <(sort cutrun_analysis/upregulated_genes.txt) | wc -l

# Real targets that are DOWN-regulated 
comm -12 <(sort geneX_real_targets.txt) <(sort cutrun_analysis/downregulated_genes.txt) | wc -l

# Decoy targets that are UP-regulated to check if it's 0 or not
comm -12 <(sort geneX_decoy_targets.txt) <(sort cutrun_analysis/upregulated_genes.txt) | wc -l

# Decoy targets that are DOWN-regulated to check if it's 0 or not
comm -12 <(sort geneX_decoy_targets.txt) <(sort cutrun_analysis/downregulated_genes.txt) | wc -l
