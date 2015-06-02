#!/bin/bash

FILES=/home/michael/batch/in/*.tif
FullProgramPath=/home/michael/Bioimage/RC_release_v1.0/RC-bactelize/B3D/RC-bactelize
OutDir=/home/michael/batch/out
ResultFilename=1-1-A_result.txt

shopt -s nullglob
for f in $FILES
do
    g=`basename "$f"`
    echo "Processing file: $f"
    $FullProgramPath "$f" -i 200 --no_fusion --no_fission --no_handle --init_mode blob_det --init_blob_min 14 --init_blob_max 24 --disc_upper -e ps --spacing 0.1802 0.1802 0.2985 --output_image "$OutDir/${g%.tif}-labeled.tif" --output_results $OutDir/$ResultFilename
done
