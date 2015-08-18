#!/bin/bash

if [ -z ${1+x} ]
then
    echo "Error!" 
    echo "Give as first parameter first part of new filename: 1- (first coverslip), 2- (second coverslip), ..." 1>&2
    echo "Give as second parameter search pattern of files to be renamed: eg.: 4-1*.ome"
    exit 1
else 

    FILES=/home/mmeuli/batch/rename/$2
    OUTDIR=/home/mmeuli/batch/in
    shopt -s nullglob
    i=1
    for f in $FILES
    do
        echo "Moving $f to $OUTDIR/$1$i.tif"
        mv "$f" "$OUTDIR/$1$i.ome"
        i=$((i+1))
    done
 
fi

