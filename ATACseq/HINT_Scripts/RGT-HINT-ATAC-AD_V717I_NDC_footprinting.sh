#!bin/bash
#This script performs TF footprinting for APP2-V717I mutation and the NDC control

#Before HINT, merge the three bam files of each mutation or control into one merged file
#First run footprinting for the V717I with the dATACpeaks
rgt-hint footprinting \
--atac-seq --paired-end \
--output-location=/home/prvaldes/scratch/ATAC/HINT/V717I \
--output-prefix=V717Ifootprint \
/home/prvaldes/scratch/ATAC/NoDuplicates/APP2-35-V717I_final_combined.bb.mapped.sorted.nd.bam \
/home/prvaldes/scratch/ATAC/BED/200923.dATACpeaks.V717I.bed \
--organism=hg38

#Now do the NDC
rgt-hint footprinting \
--atac-seq --paired-end \
--output-location=/home/prvaldes/scratch/ATAC/HINT/V717I \
--output-prefix=NDCfootprint \
/home/prvaldes/scratch/ATAC/NoDuplicates/27_4-NDC1_final_combined.bb.mapped.sorted.nd.bam \
/home/prvaldes/scratch/ATAC/BED/200923.dATACpeaks.V717I.bed \
--organism=hg38
