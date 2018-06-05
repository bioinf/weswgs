#!/bin/bash

BGLIST=$1

for i in $( cat $BGLIST )
do
	gunzip ./bedgraph/${i}.dedup.concatenated.coverage.bg.gz
	gunzip ../dedup_histograms/bedgraphs/${i}.concatenated.coverage.bg.gz
	./est_diff.sh ../dedup_histograms/bedgraphs/${i}.concatenated.coverage.bg ./bedgraph/${i}.dedup.concatenated.coverage.bg
	gzip ../dedup_histograms/bedgraphs/${i}.concatenated.coverage.bg
	gzip ./bedgraph/${i}.dedup.concatenated.coverage.bg
done
