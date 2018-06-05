#!/bin/bash

for i in histograms/*pfrag*gram ; do TAG=${i##*/} ; paste $i ../dedup_histograms/histograms/${TAG%%.dedu*}.concatenated.coverage.pfrag.histogram | awk '{ if ($2 > 0 && $1/$2 < 1) print 1-($1/$2) ; else if ($1 == 0 || $2 == 0 || $1/$2 > 1) print 0 ; else print 1 }' > ../mapble/${TAG%%.dedu*}.diff.lst ; done

cd ~/qc_proj/mandel

for i in agilent illumina roche truseq ; do paste $( grep $i NEW_VALUES_WPLS.tsv | cut -f1 | sed 's/\.dedup//' | sed 's/^/\/home\/yabarbitov\/qc_proj\/mapble\//' | perl -pe 's|\n|*diff.lst |' | perl -pe 's|$|\n|' ) | awk '{ sum=0 ; for (i=1;i<=NF;i++) sum+=$1 ; print sum/NF }' - > ${i}.diff.lst ; done
