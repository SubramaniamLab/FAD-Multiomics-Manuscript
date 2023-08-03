#!bin/bash
#This script runs motif matching for N141I vs. NDC using the HINT program

#Motif Matching of each
/Users/phoebevaldes/.local/bin/rgt-motifanalysis matching \
--organism=hg38 \
--input-files /Users/phoebevaldes/Desktop//HINT-ATAC/BED_Files/N141I_full/N141Ifootprint.bed \
--motif-dbs /Users/prvaldes/rgtdata/motifs/cisbp_hint

#Now the NDC
/Users/phoebevaldes/.local/bin/rgt-motifanalysis matching \
--organism=hg38 \
--input-files /Users/phoebevaldes/Desktop/HINT-ATAC/BED_Files/N141I_full/NDCfootprint.bed \
--motif-dbs /Users/phoebevaldes/rgtdata/motifs/cisbp_hint
