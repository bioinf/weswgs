#!/usr/bin/env python

import sys

bgf = sys.argv[1]

counter = 0
with open(bgf, 'r') as bg:
    for line in bg:
        line = line.strip()
        content = line.split('\t')
        if int(content[3]) < 10:
            counter += min(int(content[2]), int(content[6])) - max(int(content[1]), int(content[5]))

print counter
