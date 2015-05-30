#!/bin/bash

FILES=/home/michael/batch/in/*.tif
shopt -s nullglob
for f in $FILES
do
  g=`basename "$f"`
  echo $g
  echo "Processing $f file..."
  # take action on each file. $f store current file name
  /home/michael/bioimage/RC-bactelize/B/RegionCompetition "$f" -i 200 --no_fusion --no_handle --init_mode blob_det --init_blob_min 7 --init_blob_max 16 --disc_upper -e ps --spacing 0.1802 0.1802 0.2985 --output_image "/home/michael/batch/out/${g%.tif}-labeled.tif" --output_results /home/michael/batch/out/AA_result.txt -v



done
