#!/bin/bash

for i in agilent illumina roche truseq
do
	bedtools intersect -header -wa -a qc_cds_calling.GF.vcf -b ../../${i}.bed > ${i}.bait.vcf
	./numvars.py ${i}.bait.vcf Selected_Pls | grep -P "\t${i}\t" >> NUMVARS_RESTR_WGS.tsv
done
