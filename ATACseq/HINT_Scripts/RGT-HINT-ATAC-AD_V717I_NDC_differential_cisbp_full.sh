#!bin/bash

# Now run differential
# nc is number of cores
rgt-hint differential \
--organism=hg38 \
--bc \
--nc 12 \
--window-size 200 \
--standardize \
--mpbs-files=/home/prvaldes/scratch/HINT/Motif_Matching/cisBP/V717I_full/match/NDCfootprint_mpbs.bed,/home/prvaldes/scratch/HINT/Motif_Matching/cisBP/V717I_full/match/V717Ifootprint_mpbs.bed \
--reads-files=/home/prvaldes/scratch/ATAC/NoDuplicates/27_4-NDC1_final_combined.bb.mapped.sorted.nd.bam,/home/prvaldes/scratch/ATAC/NoDuplicates/APP2-35-V717I_final_combined.bb.mapped.sorted.nd.bam \
--conditions=NDC,V717I \
--output-location=/home/prvaldes/scratch/HINT/differential/cisbp/V717Ivs.NDC_full
