#!bin/bash
#This script runs FIMO
cd ~/CIS-BP
# loop for all the samples
for j in {1..1078};
	do
	FILENAME=`sed -n ''$j'p' ~/CIS-BP/cisBP_Hs_MEMEdownload_motifs.txt`
    echo $FILENAME
fimo \
--o ${FILENAME} \
--motif ${FILENAME} \
cisBP_Hs_MEMEdownload.meme \
~/genomes/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
done
