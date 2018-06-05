#!/bin/bash

# First unfiltered, second filtered
FBG=$1
SBG=$2

for tt in 0.1 0.2 0.25 0.4 0.5 0.75 0.9 0.95 0.999 1.0
do
	./count_coverage_diff.py $FBG $SBG $tt
done
