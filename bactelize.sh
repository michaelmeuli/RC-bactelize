#!/bin/bash

FILES=/home/mmeuli/batch/in/*.tif
FullProgramPath=/home/mmeuli/Bioimage/RC-bactelize/build-20150629/RC-bactelize
OutDir=/home/mmeuli/batch/out
# ResultFilename=1-1-A_result.txt

shopt -s nullglob

parallel --jobs 6 --xapply $FullProgramPath "{}" -i 200 --no_fusion --no_fission --no_handle --init_mode blob_det --init_blob_min 36 --init_blob_max 46 --disc_upper -e ps --spacing 0.0616 0.0616 0.1998 --output_image "$OutDir"/{/.}-labeled.tif --output_results "$OutDir"/{/.}-result.txt ::: $FILES


if [ $? -eq 0 ]; then
    echo OK, collecting x-x-resut.txt into 1-1-A_result.txt
    cat "$OutDir"/*.txt >> "$OutDir"/1-1-A_result.txt
else
    echo FAIL
fi





