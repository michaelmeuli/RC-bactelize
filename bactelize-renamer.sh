#!/bin/bash

EXT=h5
FILES=/media/mmeuli/WD-HD-ext4/20150903_BCG_Pasteur-Aeras_zmp1_ko_in_RAW/20160104-ijm-h5/rename/$1
OUTDIR=/media/mmeuli/WD-HD-ext4/20150903_BCG_Pasteur-Aeras_zmp1_ko_in_RAW/20160104-ijm-h5/wait

if [ -z ${1+x} ]
then
    echo "Error!" 
    echo "Give as first parameter search pattern of files to be renamed: eg.: 4-1*.ome"
    echo "Give as second parameter first part of new filename: 1- (first coverslip), 2- (second coverslip), ..." 1>&2
    exit 1

else 
    shopt -s nullglob
    i=1
    for f in $FILES
    do
        echo "Moving $f to $OUTDIR/$2$i.$EXT"
        mv "$f" "$OUTDIR/$2$i.$EXT"
        i=$((i+1))
    done
 
fi

