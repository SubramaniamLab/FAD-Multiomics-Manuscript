#!bin/bash

# Now run differential
# nc is number of cores
rgt-hint differential \
--organism=hg38 \
--bc \
--nc 12 \
--window-size 200 \
--standardize \
--mpbs-files=/home/prvaldes/scratch/HINT/Motif_Matching/cisBP/A79V_full/match/NDCfootprint_mpbs.bed,/home/prvaldes/scratch/HINT/Motif_Matching/cisBP/A79V_full/match/A79Vfootprint_mpbs.bed \
--reads-files=/home/prvaldes/scratch/ATAC/NoDuplicates/27_4-NDC1_final_combined.bb.mapped.sorted.nd.bam,/home/prvaldes/scratch/ATAC/NoDuplicates/PS1-F12424_4-A79V_final_combined.bb.mapped.sorted.nd.bam \
--conditions=NDC,A79V \
--output-location=/home/prvaldes/scratch/HINT/differential/cisbp/A79Vvs.NDC_full
