#!bin/bash
#This script removes filenames contained in the .txt file
# Loop for all the samples
while read -r filename; do
  rm "$filename"
done <empty.txt