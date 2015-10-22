#!/bin/bash

FILES=/media/mmeuli/WD-HD-ext4/20150903_BCG_Pasteur-Aeras_zmp1_ko_in_RAW/batch/in/*.ome
FullProgramPath=/home/mmeuli/Bioimage/RC-bactelize/build-20150806/RC-bactelize
OutDir=/media/mmeuli/WD-HD-ext4/20150903_BCG_Pasteur-Aeras_zmp1_ko_in_RAW/batch/out
ResultFilename=1-1-A_result.txt

shopt -s nullglob

echo -e "cNr\tiNr\tlabel\tx\ty\tz\tpSize\tpixels\tmaxDiameter\tminDiameter\troundness\tlysosome\tmacrophage" > "$OutDir"/"$ResultFilename"

parallel --jobs 6 --xapply $FullProgramPath "{}" -i 200 --no_fusion --no_fission --init_mode blob_det --init_blob_min 73 --init_blob_max 93 --disc_upper -e ps --spacing 0.0499610 0.0499610 0.1502000 --output_image "$OutDir"/{/.}-labeled.tif --output_results "$OutDir"/{/.}-result.txt ::: $FILES   #--no_handle


if [ $? -eq 0 ]; then
    echo OK, collecting x-x-resut.txt into "$ResultFilename"
    cat "$OutDir"/*.txt >> "$OutDir"/"$ResultFilename"
else
    echo FAIL
fi





