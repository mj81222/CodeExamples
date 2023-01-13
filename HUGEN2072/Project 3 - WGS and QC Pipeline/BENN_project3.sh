#!/bin/bash
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=meb301@pitt.edu
#SBATCH -t 10:00:00

set -eou pipefail

# load modules
module load fastqc/0.11.7
module load cutadapt/1.18
module load gcc/8.2.0 bwa/0.7.17
module load gcc/8.2.0 samtools/1.9
module load gatk/4.1.4.0

# fastqc
fastqc seq/hw5.fq \
       --outdir=.

# trim reads
cutadapt seq/hw5.fq \
         -m 10 \
         -q 20 \
         -f fastq \
         --quality-base=64 \
         -o reads_trimmed_proj3.fq

# fastqc trimmed
fastqc reads_trimmed_proj3.fq \
       --outdir=.

# align trimmed reads
bwa aln -t 24 \
    seq/Homo_sapiens_assembly38.fasta \
    reads_trimmed_proj3.fq \
    > reads_aligned_proj3.sai

# convert .sai to .sam file
bwa samse seq/Homo_sapiens_assembly38.fasta \
    reads_aligned_proj3.sai \
    reads_trimmed_proj3.fq \
    -f reads_aligned_proj3.sam \
    -r "@RG\tID:HW3\tLB:HW3\tSM:HW3\tPL:ILLUMINA"

# convert .sam to .bam
samtools view -bhS reads_aligned_proj3.sam > reads_aligned_proj3.bam

# mark duplicates
gatk MarkDuplicatesSpark -I reads_aligned_proj3.bam \
                         -O reads_aligned_dupsmarked_proj3.bam

# generate recalibration statistics
gatk BaseRecalibrator -I reads_aligned_dupsmarked_proj3.bam \
                      -R seq/Homo_sapiens_assembly38.fasta \
                      -O reads_BQSR_proj3.table \
                      --known-sites seq/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
                      --known-sites seq/dbsnp_146.hg38.vcf.gz

# apply recalibration statistics
gatk ApplyBQSR -R seq/Homo_sapiens_assembly38.fasta \
               -I reads_aligned_dupsmarked_proj3.bam \
               --bqsr-recal-file reads_BQSR_proj3.table \
               -O reads_aligned_dupsmarked_cleaned_proj3.bam

# generate alignment statistics
samtools flagstat reads_aligned_dupsmarked_cleaned_proj3.bam \
                  > BENN_alignment_statistics_proj3.out

# precall snps and indels
gatk HaplotypeCaller -R seq/Homo_sapiens_assembly38.fasta \
                     -I reads_aligned_dupsmarked_cleaned_proj3.bam \
                     -O genotypes_proj3.g.vcf.gz \
                     -ERC GVCF \
                     -OVI \
                     --native-pair-hmm-threads 24

# convert .g.vcf into .vcf
gatk GenotypeGVCFs -R seq/Homo_sapiens_assembly38.fasta \
                   -V genotypes_proj3.g.vcf.gz \
                   -O genotypes_proj3.vcf.gz

# create sites only .vcf
gatk MakeSitesOnlyVcf -I genotypes_proj3.vcf.gz \
                      -O sites_only_proj3.vcf.gz

# generate recalibration statistics for indels
gatk VariantRecalibrator -V sites_only_proj3.vcf.gz \
                         --trust-all-polymorphic \
                         -tranche 100.0 -tranche 99.95 -tranche 99.9 \
                          -tranche 99.5 -tranche 99.0 -tranche 97.0 \
                          -tranche 96.0 -tranche 95.0 -tranche 94.0 \
                          -tranche 93.5 -tranche 93.0 -tranche 92.0 \
                          -tranche 91.0 -tranche 90.0 \
                         -an FS -an ReadPosRankSum -an MQRankSum \
                          -an QD -an SOR -an DP \
                         -mode INDEL \
                         --max-gaussians 2 \
                         -resource:mills,known=false,training=true,truth=true,prior=12 \
                          seq/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
                         -resource:dpsnp,known=true,training=false,truth=false,prior=2 \
                          seq/dbsnp_146.hg38.vcf.gz \
                         -resource:axiomPoly,known=false,training=true,truth=false,prior=10 \
                          seq/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz \
                         -O indels_proj3.recal \
                         --tranches-file indels_proj3.tranches

# generate recalibration statistics for snps
gatk VariantRecalibrator -V sites_only_proj3.vcf.gz \
                         --trust-all-polymorphic \
                         -tranche 100.0 -tranche 99.95 -tranche 99.9 \
                          -tranche 99.8 -tranche 99.6 -tranche 99.5 \
                          -tranche 99.4 -tranche 99.3 -tranche 99.0 \
                          -tranche 98.0 -tranche 97.0 -tranche 90.0 \
                         -an QD -an MQRankSum -an ReadPosRankSum \
                          -an FS -an MQ -an SOR -an DP \
                         -mode SNP \
                         --max-gaussians 6 \
                         -resource:hapmap,known=false,training=true,truth=true,prior=15 \
                          seq/hapmap_3.3.hg38.vcf.gz \
                         -resource:omni,known=false,training=true,truth=true,prior=12 \
                          seq/1000G_omni2.5.hg38.vcf.gz \
                         -resource:1000G,known=false,training=true,truth=false,prior=10 \
                          seq/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
                         -resource:dpsnp,known=true,training=false,truth=false,prior=7 \
                          seq/dbsnp_146.hg38.vcf.gz \
                         -O snps_proj3.recal \
                         --tranches-file snps_proj3.tranches

# apply recalibration statistics of indels
gatk ApplyVQSR -V genotypes_proj3.vcf.gz \
               --recal-file indels_proj3.recal \
               --tranches-file indels_proj3.tranches \
               --truth-sensitivity-filter-level 99.7 \
               --create-output-variant-index true \
               -mode INDEL \
               -O genotypes_indelqc_proj3.vcf.gz

# apply recalibration statistics of snps
gatk ApplyVQSR -V genotypes_indelqc_proj3.vcf.gz \
               --recal-file snps_proj3.recal \
               --tranches-file snps_proj3.tranches \
               --truth-sensitivity-filter-level 99.7 \
               --create-output-variant-index true \
               -mode SNP \
               -O genotypes_indelqc_snpqc_proj3.vcf.gz

# call statistics
gatk CollectVariantCallingMetrics -I genotypes_indelqc_snpqc_proj3.vcf.gz \
                                  --DBSNP seq/dbsnp_146.hg38.vcf.gz \
                                  -O BENN_genotype_metrics_proj3
