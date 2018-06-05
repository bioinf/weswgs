#!/usr/bin/env python

import re
import sys

vcffile, plfile = sys.argv[1:]

samples = {}
samplist = []
with open(vcffile, 'r') as vf:
    for line in vf:
        line = line.strip()
        content = line.split('\t')
        if line.startswith('#'):
            if line.startswith('#CHROM'):
                for i in range(9, len(content)):
                    samples[content[i]] = [0, 0]
                    samplist.append(content[i])
            else:
                continue
        for i in range(9, len(content)):
            sampl = content[i].split(':')
            if (sampl[0] == '0/1' or sampl[0] == '1/1') and len(content[3]) != len(content[4]):
               samples[samplist[i - 9]][0] += 1
               if 'lowGQ' in content[i]:
                   samples[samplist[i - 9]][1] += 1

platforms =  {}
with open(plfile, 'r') as pf:
    for line in pf:
        content = line.strip().split('\t')
        platforms[content[0]] = content[1]

for jj in platforms:
    if jj not in samples:
        continue
    if 'sample_524' in jj:
        continue
    ars = samples[jj]
    pl = platforms[jj]
    print str(ars[0]) + '\t' + pl + '\t' + 'ALL'
    print str(ars[1]) + '\t' + pl + '\t' + 'LOWGQ'
