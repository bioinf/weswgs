#!/bin/bash

bedtools multicov -bams ../cds_calling/bams/dedupped/*.bam -bed ./CLV_Deciles_normal.bed > sites_coverage.bed
for i in $( ls ../cds_calling/bams/dedupped/*.bam | grep -oP 'sample_\K[^\.]+' ) ; do grep -P "[S_]${i}[_\.]" ../final.list ; done | perl -pe 's/.*_//g' > pls.list
for i in `seq 1 191678` ; do cat pls.list ; done > platforms.col
awk '{OFS="\t"} { for (i=5;i<=NF;i++) print $4,$i }' sites_coverages.bed > coverages.col
paste coverages.col platforms.col > ALL_NEW_GCBIAS.hist
