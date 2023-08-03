#!bin/bash
#This script merges three BAM files of each N141I mutation into one merged BAM file.
#For ATAC Samples
module load samtools
cd /home/prvaldes/scratch/ATAC/NoDuplicates
#Perform the merge command
samtools merge PS2-F12448_9-A141I_final_combined.bb.mapped.sorted.nd.bam PS2-F12448_9-A141I-1_S13_combined.bb.mapped.sorted.nd.bam PS2-F12448_9-A141I-2_S14_combined.bb.mapped.sorted.nd.bam PS2-F12448_9-A141I-3_S15_combined.bb.mapped.sorted.nd.bam