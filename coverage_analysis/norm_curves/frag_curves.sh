#!/bin/bash


for i in *.pfrag*
do
	MCOV=$( awk '{ sum+=$1 } END { print sum/NR }' $i )
        PL=$( grep ${i%%.concat*} ../../final.list | grep -oP '_\K[arit].*' )
        ./make_pfrag_curve.py $i $MCOV $PL
done
