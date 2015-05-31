#!/bin/bash

if [ -z ${1+x} ]
then
    echo "Error! Give first part of new filename as parameter (1A, 1B, ...)" 1>&2
    exit 1
else 

    FILES=/home/michael/batch/tif-orig/*.tif
    OUTDIR=/home/michael/batch/in
    shopt -s nullglob
    i=1
    for f in $FILES
    do
        echo "Moving $f to $OUTDIR/$1$i"
        mv "$f" "$OUTDIR/$1$i"
        i=$((i+1))
    done
 
fi

