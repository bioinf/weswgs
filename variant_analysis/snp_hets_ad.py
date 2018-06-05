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
                    samples[content[i]] = []
                    samplist.append(content[i])
            else:
                continue
        if len(content[3]) != len(content[4]):
            continue
#        print samplist
#        print samples
        formatstr = content[8].split(':')
        if 'AD' not in formatstr:
            continue
        ad_index = formatstr.index('AD')
        for i in range(9, len(content)):
            sampl = content[i].split(':')
            if 'lowGQ' in content[i]:
                continue
            if ad_index >= len(sampl) or '.' in sampl[ad_index] or sampl[0] != '0/1':
                continue
            ref_allele = int(sampl[ad_index].split(',')[0])
            alt_allele = int(sampl[ad_index].split(',')[1])
            if ref_allele == 0 or alt_allele == 0:
                continue
            allele_ratio = alt_allele/float((ref_allele + alt_allele))
            samples[samplist[i - 9]].append(allele_ratio)

platforms =  {}
with open(plfile, 'r') as pf:
    for line in pf:
        content = line.strip().split('\t')
        platforms[content[0]] = content[1]

for jj in platforms:
    ars = samples[jj]
    pl = platforms[jj]
    for ar in ars:
        print str(ar) + '\t' + pl    
