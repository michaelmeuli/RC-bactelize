#!/bin/bash

FILES=/home/michael/batch/in/*.tif
FullProgramPath=/home/michael/bioimage/RC-bactelize/B/RC-bactelize
OutDir=/home/michael/batch/out
ResultFilename=AA_result.txt

shopt -s nullglob
for f in $FILES
do
    g=`basename "$f"`
    echo "Processing file: $f"
    $FullProgramPath "$f" -i 200 --no_fusion --no_handle --init_mode blob_det --init_blob_min 7 --init_blob_max 16 --disc_upper -e ps --spacing 0.1802 0.1802 0.2985 --output_image "$OutDir/${g%.tif}-labeled.tif" --output_results $OutDir/$ResultFilename --verbose
done
