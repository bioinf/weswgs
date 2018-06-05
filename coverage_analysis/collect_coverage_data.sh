#!/bin/bash

BAM=$1
INTERVALS=$2

echo "Processing $BAM"
bedtools genomecov -ibam $BAM -bga | bedtools intersect -wo -a - -b $INTERVALS | sort -n -k6,6 -n -k7,7 -n -k2,2 -n -k3,3 - > ${BAM%%.dedup.bam}.concatenated.coverage.bg 
wait

./bg_to_pfcov.py ${BAM%%.dedup.bam}.concatenated.coverage.bg
MCOV=$( awk '{ sum+=$1 } END { print sum/NR }' ${BAM%%.dedup.bam}.concatenated.coverage.pfrag.histogram )
LOWCOV=$( ./length_10x.py ${BAM%%.dedup.bam}.concatenated.coverage.bg )

echo ${BAM%%.dedup*} $MCOV $LOWCOV
