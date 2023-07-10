# quick script to take output of verifySplicing and convert to a list of one isoform per line. This resulting list is then used in the consolidate_isos.py script.

import sys

with open(sys.argv[1]) as infile:
    for line in infile:
        newline = line.strip().split()
        isos = newline[-1].split(",")
        for iso in isos:
            print(iso)
