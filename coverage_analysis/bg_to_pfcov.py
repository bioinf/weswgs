#!/usr/bin/env python

import numpy
import sys

bgfile = sys.argv[1]

outprefix = bgfile.rstrip('.bg')

with open(bgfile, 'r') as bgf, open(outprefix + '.pbase.histogram', 'w') as pbh, open(outprefix + '.pfrag.histogram', 'w') as pfh:
    this_bait = (0, 0)
    for line in bgf:
        content = line.split('\t')
        if (int(content[5]), int(content[6])) != this_bait:
            if this_bait != (0, 0):
                mean_cov = numpy.mean(frag_cov)
                pfh.write(str(mean_cov) + '\n')
                if mean_cov >= 20 and len(frag_cov) == 100:
                    norm_cov = [str((x * 1.0)/mean_cov) for x in frag_cov]
                    pbh.write(str(this_bait[1] - this_bait[0]) + '\t' + '\t'.join(norm_cov) + '\n')
            this_bait = (int(content[5]), int(content[6]))
            breaks = numpy.linspace(this_bait[0], this_bait[1], 100)
            counter = 0
            frag_cov = []
        while counter <= 99:
            if int(content[1]) <= breaks[counter] <= int(content[2]):
                frag_cov.append(int(content[3]))
                counter += 1
            elif breaks[counter] > int(content[2]):
                break
#            else:
#                print 'Looping over bullshit', counter, content[1], content[2], ' on bait ', this_bait
