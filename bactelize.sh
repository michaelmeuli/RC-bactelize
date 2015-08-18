#!/bin/bash

FILES=/home/mmeuli/batch/in/*.ome
FullProgramPath=/home/mmeuli/Bioimage/RC-bactelize/build-20150806/RC-bactelize
OutDir=/home/mmeuli/batch/out
ResultFilename=1-1-A_result.txt

shopt -s nullglob

echo -e "coverslipNr\timageNr\tlabel\tx\ty\tz\tphysicalSize\tnumberOfPixels\tmaxDiameter\tminDiameter\troundness\tmean_lysosome\tmean_macrophage" > "$OutDir"/"$ResultFilename"

parallel --jobs 6 --xapply $FullProgramPath "{}" -i 200 --no_fusion --no_fission --init_mode blob_det --init_blob_min 36 --init_blob_max 46 --disc_upper -e ps --spacing 0.0616 0.0616 0.1998 --output_image "$OutDir"/{/.}-labeled.tif --output_results "$OutDir"/{/.}-result.txt ::: $FILES   #--no_handle


if [ $? -eq 0 ]; then
    echo OK, collecting x-x-resut.txt into "$ResultFilename"
    cat "$OutDir"/*.txt >> "$OutDir"/"$ResultFilename"
else
    echo FAIL
fi





