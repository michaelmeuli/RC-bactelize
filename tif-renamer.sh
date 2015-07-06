#!/bin/bash

if [ -z ${1+x} ]
then
    echo "Error! Give first part of new filename as parameter: 1- (first coverslip), 2- (second coverslip), ..." 1>&2
    exit 1
else 

    FILES=/home/mmeuli/batch/tif-orig/*.tif
    OUTDIR=/home/mmeuli/batch/in
    shopt -s nullglob
    i=1
    for f in $FILES
    do
        echo "Moving $f to $OUTDIR/$1$i.tif"
        mv "$f" "$OUTDIR/$1$i.tif"
        i=$((i+1))
    done
 
fi

