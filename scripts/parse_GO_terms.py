# This script is used to parse (collapse) Gene Ontology association and to convert to a format used by goatools python package.

import pandas as pd #use pandas to convert dictionary to dataframe at end

#file = open('/data3/peter/sunflower/ref_genome_Ha412HO/HAN412_vs_ATH_GO_GOSLIM.uniq.txt', 'r')
file = open('/data3/peter/sunflower/ref_genome_Araport11/tmp', 'r')
parsed_GO={}
prev_gene='AT1G01020' #set prev gene as the first gene in list to start
#prev_GO_term='GO:0003674'
for current_line in file:
    
    current_gene = current_line.rstrip().split('\t')[0]
    current_GO_term = current_line.rstrip().split('\t')[1]
    if current_gene not in parsed_GO: #if it's a new gene
        parsed_GO[current_gene]=[] #then make a new list
        parsed_GO[current_gene].append(current_GO_term) #and add the GO term
    else: parsed_GO[current_gene].append(current_GO_term) #if same gene as previous line, just add the GO term
    
    prev_gene=current_gene #reset prev gene

#print(parsed_GO)

# print as a tab separated table
for k,v in parsed_GO.items():
    print(k+"\t",v)

