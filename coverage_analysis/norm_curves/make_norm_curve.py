#!/usr/bin/env python

import numpy as np
import re
import sys
import os

bgfi, mcov, plf = sys.argv[1:]

mcov = float(mcov)

normcovs = {}
points = np.linspace(0.00, 3.0, 301)
for jj in points:
    normcovs["%.2f" % jj] = 0

with open(bgfi, 'r') as bgf:
    for line in bgf:
        content = line.split('\t')
        leftborder = np.max([int(content[5]), int(content[1])])
        rightborder = np.min([int(content[6]), int(content[2])])
        lencov = rightborder - leftborder
        cov = float(content[3])
        ncov = cov/mcov
        if ncov < 3.0:
            normcovs["%.2f" % ncov] += lencov
        else:
            normcovs['3.00'] += lencov

noint = sum(normcovs.values())
remcov = noint
for i in points:
    print "%.2f" % i + '\t' + str((remcov - normcovs["%.2f" % i])/(noint * 1.0)) + '\t' + plf
    remcov -= normcovs["%.2f" % i]
