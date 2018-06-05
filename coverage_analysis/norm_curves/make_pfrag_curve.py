#!/usr/bin/env python

import sys
import numpy as np

vcfi, mcov, plf = sys.argv[1:]
mcov = float(mcov)

normcovs = {}
points = np.linspace(0.00, 3.0, 301)
for jj in points:
    normcovs["%.2f" % jj] = 0

with open(vcfi, 'r') as vf:
    for line in vf:
        cov = float(line.strip())
        ncov = cov/mcov
        if ncov < 3.0:
            normcovs["%.2f" % ncov] += 1
        else:
            normcovs['3.00'] += 1

#print sum(normcovs.values())
noint = sum(normcovs.values())
remcov = noint
for i in points:
    print "%.2f" % i + '\t' + str((remcov - normcovs["%.2f" % i])/(noint * 1.0)) + '\t' + plf
    remcov -= normcovs["%.2f" % i]
