#!bin/bash
#This script performs TF footprinting for PSEN1-A79V mutation and the NDC control

#Before HINT, merge the three bam files of each mutation or control into one merged file
#First run footprinting for the A79V with the dATACpeaks
rgt-hint footprinting \
--atac-seq --paired-end \
--output-location=/home/prvaldes/scratch/ATAC/HINT/A79V \
--output-prefix=A79Vfootprint \
/home/prvaldes/scratch/ATAC/NoDuplicates/PS1-F12424_4-A79V_final_combined.bb.mapped.sorted.nd.bam \
/home/prvaldes/scratch/ATAC/BED/200923.dATACpeaks.A79V.bed \
--organism=hg38

#Now do the NDC
rgt-hint footprinting \
--atac-seq --paired-end \
--output-location=/home/prvaldes/scratch/ATAC/HINT/A79V \
--output-prefix=NDCfootprint \
/home/prvaldes/scratch/ATAC/NoDuplicates/27_4-NDC1_final_combined.bb.mapped.sorted.nd.bam \
/home/prvaldes/scratch/ATAC/BED/200923.dATACpeaks.A79V.bed \
--organism=hg38
