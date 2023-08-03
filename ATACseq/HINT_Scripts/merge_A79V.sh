#!bin/bash
#This script merges three BAM files of each A79V mutation into one merged BAM file.
#For ATAC Samples
module load samtools
cd /home/prvaldes/scratch/ATAC/NoDuplicates
#Perform the merge command
samtools merge PS1-F12424_4-A79V_final_combined.bb.mapped.sorted.nd.bam PS1-F12424_4-A79V-1_S22_combined.bb.mapped.sorted.nd.bam PS1-F12424_4-A79V-2_S23_combined.bb.mapped.sorted.nd.bam PS1-F12424_4-A79V-3_S24_combined.bb.mapped.sorted.nd.bam