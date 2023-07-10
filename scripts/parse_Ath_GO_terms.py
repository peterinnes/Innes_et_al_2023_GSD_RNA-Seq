# This script is used to parse (collapse) Gene Ontology association and to convert to a format used by goatools python package.
# ATH_GO_GOSLIM.txt downloaded from arabidopsis.org
# input file 'tmp' generated with:
# tail -n +5 ATH_GO_GOSLIM.txt | cut -f1,6 | sort -k1,1 -k2,2 | uniq > tmp
# also ended up parsing the output of this python script to produce the final 'id2gos' format:
# python3.6 ~/gsd_RNA-seq/scripts/parse_GO_terms.py | sed "s/'//g" | sed 's/,/;/g' | tr -d '\ []' > ~/gsd_RNA-seq/data/ref_genome_Araport11/ATH_GO_GOSLIM.id2gos.txt

import pandas as pd #use pandas to convert dictionary to dataframe at end

file = open('~/gsd_RNA-seq/data/ref_genome_Araport11/tmp', 'r')
parsed_GO={}
prev_gene=''
for current_line in file:
    
    current_gene = current_line.rstrip().split('\t')[0]
    current_GO_term = current_line.rstrip().split('\t')[1]
    if current_gene not in parsed_GO: #if it's a new gene
        parsed_GO[current_gene]=[] #then make a new list
        parsed_GO[current_gene].append(current_GO_term) #and add the GO term
    else: parsed_GO[current_gene].append(current_GO_term) #if same gene as previous line, just add the GO term
    
    prev_gene=current_gene #reset prev gene

# print as a tab separated table
for k,v in parsed_GO.items():
    print(k+"\t",v)

