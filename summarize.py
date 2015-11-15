#!/usr/bin/env python

import sys
import os

results = []

result_files = [x for x in os.listdir(sys.argv[1]) if x.endswith(".tsv")]

for f in result_files:
    for line in open(f, 'r'):
        if line.startswith(">"):
            line = line[1:]
            line = line.split()
            results.append((line[0], line[1]))

for result in results:
    print >> sys.stdout, "{0}\t{1:.4f}".format(result[0], float(result[1]))

