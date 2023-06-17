#!/usr/bin/env python

import sys, os, re
from collections import defaultdict



usage = "\n\n\tusage: {} chims.described > FI.list\n\n".format(sys.argv[0])

if len(sys.argv) < 2:
    exit(usage)

chims_described_file = sys.argv[1]

fusion_counter = defaultdict(int)

with open(chims_described_file) as fh:
    for line in fh:
        line = line.rstrip()
        vals = line.split(";")
        fusion_name = vals.pop()
        fusion_counter[fusion_name] += 1

fusions = fusion_counter.keys()
fusions = sorted(fusions, key=lambda x: fusion_counter[x], reverse=True)

for fusion in fusions:
    print("\t".join([fusion, str(fusion_counter[fusion])]))


sys.exit(0)

