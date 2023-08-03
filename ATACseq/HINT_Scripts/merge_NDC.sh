#!bin/bash
#This script merges three BAM files of each control into one merged BAM file.
#For ATAC Samples
module load samtools
cd /home/prvaldes/scratch/ATAC/NoDuplicates
#Perform the merge command
samtools merge 27_4-NDC1_final_combined.bb.mapped.sorted.nd.bam 27_4-NDC1-1_S1_combined.bb.mapped.sorted.nd.bam 27_4-NDC1-2_S2_combined.bb.mapped.sorted.nd.bam 27_4-NDC1-3_S3_combined.bb.mapped.sorted.nd.bam
