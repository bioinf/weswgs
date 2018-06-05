#!/usr/bin/env python

import sys
import re

firstbg, secondbg, threshold = sys.argv[1:]
threshold = float(threshold)

with open(firstbg, 'r') as fbg:
    coverages_f = []
    for line in fbg:
        content = line.strip().split('\t')
        lborder = max(int(content[5]), int(content[1]))
        rborder = min(int(content[6]), int(content[2]))
        for kk in range(rborder - lborder):
            coverages_f.append(int(content[3]))

with open(secondbg, 'r') as sbg:
    coverages_s = []
    for line in sbg:
        content = line.strip().split('\t')
        lborder = max(int(content[5]), int(content[1]))
        rborder = min(int(content[6]), int(content[2]))
        for kk in range(rborder - lborder):
            coverages_s.append(int(content[3]))

diffcounter = 0
for kk in range(len(coverages_s)):
    if coverages_s[kk] <= ((1 - threshold) * coverages_f[kk]):
        diffcounter += 1

print re.findall('wes_\d+.sample_[^\.]+', firstbg)[0] + '\t' + str(threshold) + '\t' + str(diffcounter)
