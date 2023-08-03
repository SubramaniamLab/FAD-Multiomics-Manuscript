#!bin/bash
#This script replaces underscores with . in .bed files
# loop for all the samples and iterate over the filenames
for f in *_*.*_*.bed; do i="${f%.bed}"; mv -i -- "$f" "${i//_/.}.bed"; done
