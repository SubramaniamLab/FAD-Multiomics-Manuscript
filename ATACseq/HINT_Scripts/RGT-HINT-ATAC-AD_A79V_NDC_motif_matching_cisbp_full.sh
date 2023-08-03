#!bin/bash
#This script runs motif matching for A79V vs. NDC using the HINT program

#Motif Matching of each
/Users/phoebevaldes/.local/bin/rgt-motifanalysis matching \
--organism=hg38 \
--input-files /Users/phoebevaldes/Desktop/HINT-ATAC/BED_Files/A79V_full/A79Vfootprint.bed \
--motif-dbs /Users/prvaldes/rgtdata/motifs/cisbp_hint

#Now the NDC
/Users/phoebevaldes/.local/bin/rgt-motifanalysis matching \
--organism=hg38 \
--input-files /Users/phoebevaldes/Desktop/HINT-ATAC/BED_Files/A79V_full/NDCfootprint.bed \
--motif-dbs /Users/phoebevaldes/rgtdata/motifs/cisbp_hint
