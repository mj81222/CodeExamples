#!/bin/bash
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=meb301@pitt.edu
#SBATCH -t 10:00:00

set -eou pipefail

# load packages
module load gcc/8.2.0 samtools/1.9
module load samtools/1.12

# set directory
cd ~/2072/project4

# extract chr22 and convert to .bam format
samtools view -T /bgfs/rminster/hugen2072-2021s/data/P4/GRCh38_full_analysis_set_plus_decoy_hla.fa \
              -bS \
              -o ./NA12778.chr22.bam  \
              -h /bgfs/rminster/hugen2072-2021s/data/P4/NA12778.final.cram chr22

# index the file
samtools index ./NA12778.chr22.bam

# set up manta pipeline
python /ihome/crc/install/manta/manta-1.6.0.centos6_x86_64/bin/configManta.py \
       --bam ./NA12778.chr22.bam \
       --referenceFasta /bgfs/rminster/hugen2072-2021s/data/P4/GRCh38_full_analysis_set_plus_decoy_hla.fa \
       --runDir manta_test

# run manta
## change to manta directory
cd ./manta_test

## run the script
python runWorkflow.py

# results
cd results/variants

# extract deletions and duplications and convert to .bed
zcat diploidSV.vcf.gz \
  | awk '{if ($5~"DEL" || $5~"DUP") print}'\
  | awk -F"[\t;]" '{print $1,$2,$8,$5}' OFS="\t" \
  | sed 's/END=//g'> cnv.bed

# extract larger than 1kb
awk '{if ( $3-$2>=1000) print}' cnv.bed > BENN_gt1kb.cnv.bed

# change directory
cd ~/2072/project4

# convert to .bam file and index
samtools view -h -o BENN_NA12778.chr22.cram NA12778.chr22.bam
samtools index BENN_NA12778.chr22.cram
