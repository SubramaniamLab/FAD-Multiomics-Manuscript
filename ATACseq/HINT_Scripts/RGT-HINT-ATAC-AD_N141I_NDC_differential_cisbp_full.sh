#!bin/bash

# Now run differential
# nc is number of cores
rgt-hint differential \
--organism=hg38 \
--bc \
--nc 12 \
--window-size 200 \
--standardize \
--mpbs-files=/home/prvaldes/scratch/HINT/Motif_Matching/cisBP/N141I_full/match/NDCfootprint_mpbs.bed,/home/prvaldes/scratch/HINT/Motif_Matching/cisBP/N141I_full/match/N141Ifootprint_mpbs.bed \
--reads-files=/home/prvaldes/scratch/ATAC/NoDuplicates/27_4-NDC1_final_combined.bb.mapped.sorted.nd.bam,/home/prvaldes/scratch/ATAC/NoDuplicates/PS2-F12448_9-A141I_final_combined.bb.mapped.sorted.nd.bam \
--conditions=NDC,N141I \
--output-location=/home/prvaldes/scratch/HINT/differential/cisbp/N141Ivs.NDC_full
