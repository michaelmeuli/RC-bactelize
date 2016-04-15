#!/bin/bash

OutDir=/media/mmeuli/WD-HD-ext4/20160408_BCG_Pasteur-Aeras_in_THP-1/data-deconvolved/out
ResultFilename=1-1-A_result.txt

shopt -s nullglob

echo -e "cNr\tiNr\tlabel\tx\ty\tz\tpSize\tpixels\tmaxDiameter\tminDiameter\troundness\tlysosome\tmacrophage" > "$OutDir"/"$ResultFilename"

echo OK, collecting x-x-resut.txt into "$ResultFilename"
cat "$OutDir"/*-result.txt >> "$OutDir"/"$ResultFilename"






