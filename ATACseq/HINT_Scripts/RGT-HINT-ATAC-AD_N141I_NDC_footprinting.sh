#!bin/bash
#This script performs TF footprinting for PSEN2-N141I mutation and the NDC control

#Before HINT, merge the three bam files of each mutation or control into one merged file
#First run footprinting for the N141I with the dATACpeaks
rgt-hint footprinting \
--atac-seq --paired-end \
--output-location=/home/prvaldes/scratch/ATAC/HINT/N141I \
--output-prefix=N141Ifootprint \
/home/prvaldes/scratch/ATAC/NoDuplicates/PS2-F12448_9-A141I_final_combined.bb.mapped.sorted.nd.bam \
/home/prvaldes/scratch/ATAC/BED/200923.dATACpeaks.A141I.bed \
--organism=hg38

#Now do the NDC
rgt-hint footprinting \
--atac-seq --paired-end \
--output-location=/home/prvaldes/scratch/ATAC/HINT/N141I \
--output-prefix=NDCfootprint \
/home/prvaldes/scratch/ATAC/NoDuplicates/27_4-NDC1_final_combined.bb.mapped.sorted.nd.bam \
/home/prvaldes/scratch/ATAC/BED/200923.dATACpeaks.A141I.bed \
--organism=hg38
