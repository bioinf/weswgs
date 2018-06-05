#!/bin/bash


for i in *.bg.gz
do
	gunzip $i
	MCOV=$( awk '{ sum+=$1 } END { print sum/NR }' ../histograms/${i%%.bg.gz}.pfrag.histogram )
        PL=$( grep ${i%%.concat*} ../../mandel/NEW_VALUES_WPLS.tsv | cut -f 4 )
        ./make_norm_curve.py ${i%%.gz} $MCOV $PL
	gzip ${i%%.gz}
done
