#!bin/bash
#This script merges three BAM files of each V717I mutation into one merged BAM file.
#For ATAC Samples
module load samtools
cd /home/prvaldes/scratch/ATAC/NoDuplicates
#Perform the merge command
samtools merge APP2-35-V717I_final_combined.bb.mapped.sorted.nd.bam APP2-35-V717I-1_S10_combined.bb.mapped.sorted.nd.bam APP2-35-V717I-2_S11_combined.bb.mapped.sorted.nd.bam APP2-35-V717I-3_S12_combined.bb.mapped.sorted.nd.bam