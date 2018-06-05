#!/bin/bash

PL=$1

for i in $( cat $PL )
do
	awk '{ for(i=1;i<=NF;i++) total[i]+=$i ; } END { for(i=1;i<=NF;i++) printf "%f ",total[i]/NR ;}' ./raw/${i}.concatenated.coverage.pbase.histogram | perl -pe 's|$|\n|'
done
