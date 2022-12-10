#Outputs list of the longest isoform of each gene, along with its length.
#example command: python get_longest_iso.py transcriptome.fasta > longest_iso_list.txt

import sys, numpy
from operator import itemgetter

#### read in transcriptome and find longest iso ####
lengthDictionary = {}
with open(sys.argv[1], "r") as infile:
    for line in infile:
        if line[0] == ">":    
            newline = line.strip().split()
            gene = "_".join(newline[0].split(">")[1].split("_")[0:4])
            iso = newline[0].split(">")[1].split("_")[4]
        else:
            length = len(line.strip())
            if gene not in lengthDictionary:
                lengthDictionary[gene] = [iso, length]
            else:
                if length > lengthDictionary[gene][1]:
                    lengthDictionary[gene] = [iso, length]

for gene in lengthDictionary.keys():
    print gene+'_'+lengthDictionary[gene][0],'\t',lengthDictionary[gene][1]
