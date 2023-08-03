#!bin/bash
#This script converts bed name
# loop for all the samples
for j in {1..988};
	do
	FILENAME=`sed -n ''$j'p' CISBP_human_TFnames.txt`
    echo $FILENAME
mv ${FILENAME}.bed ${FILENAME}_TFBS.bed
done
