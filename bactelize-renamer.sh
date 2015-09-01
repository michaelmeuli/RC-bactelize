#!/bin/bash

if [ -z ${1+x} ]
then
    echo "Error!" 
    echo "Give as first parameter search pattern of files to be renamed: eg.: 4-1*.ome"
    echo "Give as second parameter first part of new filename: 1- (first coverslip), 2- (second coverslip), ..." 1>&2
    exit 1
else 

    EXT=ome
    FILES=/home/mmeuli/batch/rename/$1
    OUTDIR=/home/mmeuli/batch/in
    shopt -s nullglob
    i=1
    for f in $FILES
    do
        echo "Moving $f to $OUTDIR/$2$i.$EXT"
        mv "$f" "$OUTDIR/$2$i.$Ext"
        i=$((i+1))
    done
 
fi

