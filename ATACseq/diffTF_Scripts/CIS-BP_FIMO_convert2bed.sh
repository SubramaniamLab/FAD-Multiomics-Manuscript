#!bin/bash
#This script converts gff3 to bed
cd ~/CIS-BP/CIS-BP_MEME_gff3
# loop for all the samples
for j in {1..988};
	do
	FILENAME=`sed -n ''$j'p' ~/CIS-BP/CISBP_human_TFnames.txt`
    echo $FILENAME
cd ~/CIS-BP/CIS-BP_MEME_gff3/${FILENAME}
gff2bed < fimo.gff > ${FILENAME}_full.bed
awk '{ gsub("qvalue= ","qvalue="); print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10 }' ${FILENAME}_full.bed > ${FILENAME}_corrected.bed
awk '{split($10,a,"sequence=");print $1,$2,$3,a[2],$5,$6}' OFS='\t' ${FILENAME}_corrected.bed > ${FILENAME}_int.bed
tr -d ';' < ${FILENAME}_int.bed > ~/CIS-BP/CIS-BP_MEME_TFBS_human/${FILENAME}.bed
rm ${FILENAME}_corrected.bed
rm ${FILENAME}_int.bed
done
