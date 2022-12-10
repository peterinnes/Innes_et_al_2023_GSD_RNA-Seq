################
#### PARENTS_DIFF_V2.PY this program identifies alternatively-spliced ####
#### genes with significantly different isoform composition between ####
#### ecotypes. Significance is assesses either by t-test or MANOVA ####
#### if there are more that two isoforms. Isoform composition is ####
#### transformed by isometric log ratio (ILR) prior to significance test. ####
#### example command: python3 parents_diff_v2 good_isoforms.txt consolidated_isoforms.txt isoform_quant_list.txt ####
#### isoform quant list is a list of file paths to raw isoform count data of each sample ####

import sys, os
from operator import itemgetter
import numpy as np
import pandas as pd
from skbio.stats.composition import ilr
from statsmodels.multivariate.manova import MANOVA
from statsmodels.multivariate.pca import PCA
from scipy import stats

#trans_path = sys.argv[1] #commenting this out because we aren't looking at transgressive isoforms.
goodisos_path = sys.argv[1]
consolidation_path = sys.argv[2]
parent_list = sys.argv[3]



##### read in list of transgressive isoforms #####
#list_of_trans_isos = {}
#list_of_trans_genes = {}
#with open(trans_path) as infile:
#    for line in infile:
#        newline = line.strip().split()
#        myiso = newline[0]
#        mygene = "_".join(myiso.split("_")[0:4])
#        print(myiso, mygene)
#        list_of_trans_isos[myiso] = 0
#        if mygene not in list_of_trans_genes:
#            list_of_trans_genes[mygene] = []
#        # unindent
#        list_of_trans_genes[mygene].append( myiso )
#







##### read in list of good isoforms. 'new_gene' refers to our new notation using decimals to denote new gene groups (those that were idenitified as a single gene by trinity, but aligned to different regions of the genome). 'old_gene' is the original gene notation #####   
good_isos = {}
with open(goodisos_path) as infile:
    for line in infile:
        new_gene, isos = line.strip().split()
        isos = isos.split(",")
        for iso in isos:
            old_gene = new_gene.split(".")[0]
            if old_gene not in good_isos: #some are repeated
                good_isos[old_gene] = {}
            good_isos[old_gene][new_gene] = []
            for iso in isos:
                good_isos[old_gene][new_gene].append(iso)
        # old code
        #gene = gene.split(".")[0]
        #good_isos[gene] = []
        #for iso in isos:
            #good_isos[gene].append(iso)

# Commenting this out bevause I don't think it's necessary without a list of trans_isos
#        for iso in isos:
#            if iso in list_of_trans_isos:
#                old_gene = new_gene.split(".")[0]
#                #print(old_gene)
#                if old_gene not in good_isos: # some are repeated
#                    good_isos[old_gene] = {}
#                good_isos[old_gene][new_gene] = [] # adding in the new gene
#                for iso in isos:
#                    good_isos[old_gene][new_gene].append(iso)
#


                    

##### read in list of consolidated alleles #####
consolidation = {}
with open(consolidation_path) as infile:
    for line in infile:
        iso_ID, others = line.strip().split()
        consolidation["_".join(iso_ID.split("_")[0:4])] = []
        for other in others.split("|"):
            for indiv in other.split(","):
                consolidation[indiv] = other.split(",")

            


##### read in parent expression data. "parents" in our case are just all the dune/non-dune samples#####
sys.stderr.write("reading in individuals' tpm data\n")
parents = {} # dictionary of parents (samples), genes within parents, isoforms within genes
with open(parent_list) as list_file:
    for line in list_file:
        parent_path = line.strip()
        parent_id = parent_path.split("/")[-1].split(".")[-3]
        #print("getting tpm data for", parent_id)
        parents[parent_id] = {}
        with open(line.strip()) as parent_file:
            parent_file.readline() # header
            for line in parent_file:
                newline = line.strip().split() 
                iso = newline[0]
                old_gene = "_".join(iso.split("_")[0:4])
                if old_gene in good_isos:
                    tpm = float(newline[5]) #changed from newline[3] in Chris' original code to work with slightly different file format
                    if old_gene not in parents[parent_id]:
                        parents[parent_id][old_gene] = {}
                    # unindent
                    parents[parent_id][old_gene][iso] = tpm




##### now the business #####
sys.stderr.write("calculating differential splicing...")
for old_gene in good_isos:
    for new_gene in good_isos[old_gene]:
        #print("CURRENT GENE is", new_gene)
        #### Calculate total expression of current gene for each parent ####
        parent_ids = parents.keys()
        totals = [0] * len(parent_ids) #list containing expression level of current gene in each parent/sample
        parent_counter=0
        iso_list = good_isos[old_gene][new_gene]
        for parent in parents:
            #print("totalling expression of current gene for", parent)
            for iso in iso_list:
                #print(iso, parents[parent][gene][iso])
                totals[parent_counter] += parents[parent][old_gene][iso] + 0.000001 #small value added to avoid zeros. 'totals' is a list so we need to use a counter to index, as the loop is otherwise through the 'parents' dictionary. was getting an error with original code so did this instead. there is probably a less janky fix.
            parent_counter += 1
     

        ##### calculate proportions #####
        #print("finding alleles of the same splice form and calculating expression proportions...")
        isoform_clusters = []
        for iso in iso_list:
            if iso in consolidation:
                if consolidation[iso] not in isoform_clusters:
                    isoform_clusters.append(consolidation[iso])
            else:
                isoform_clusters.append([iso])
        #print(isoform_clusters)
        props = [ [0]*len(isoform_clusters) for i in range(len(parent_ids))]
        
        parent_counter=0
        for parent in parents:
            indiv_props = []
            cluster_totals = [0]*len(isoform_clusters)
            for cluster in range(len(isoform_clusters)):
                for allele in isoform_clusters[cluster]: #sum up the alleles of the same splice form
                    cluster_totals[cluster] += parents[parent][old_gene][allele] + 0.000001
                # unindent
                props[parent_counter][cluster] =  cluster_totals[cluster] / totals[parent_counter]  #then divide by the total gene expression for that parent
            parent_counter += 1

        count_isos_in_parents = 0 # real quick
        for cluster in range(len(isoform_clusters)):
            found = False
            for parent in range(len(parent_ids)):
                if props[parent][cluster] > 1e-05: #setting threshold for cluster to be counted or not.
                    found = True
            if found == True:
                count_isos_in_parents += 1
        #print("allele clusters found in", old_gene, "across all samples:", count_isos_in_parents)
        #print("isoform proportions:", props)      


        ##### Isometric Log Ratio  transform #####
        ilrs = [None]*len(parent_ids)
        for parent in range(len(parent_ids)):
            if 1 not in props[parent]: #hacky fix to get around proportions of 1. 
                ilrs[parent] = ilr(props[parent]) #this breaks if 'props' contains proportions of 1. but I don't think we should consider proportions of 1, bc then splicing isn't occuring? 
        ilrs = np.array(ilrs)
        #print("isometric log raio transform of proportions:", ilrs)


        ##### test #####
        if None not in ilrs: #avoiding proportions of 1, again.
            if len(ilrs.shape) > 1: #if more than one ILR column, use manova
                pca_model = PCA(ilrs)
                threshold = 0.99
                if pca_model.rsquare[1] < threshold: #filter for collinearity 
                    pops = np.array([1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0]) # should encode for population of each parent/sample ['dune-10', 'non_dune-16']
                    manova = MANOVA(endog=ilrs, exog=pops)
                    #print(manova.mv_test())
                    manova_pval = manova.mv_test().results['x0']['stat'].values[0,4]
                    print(new_gene, manova_pval)
                    #if manova_pval < .05:
                    #    print(new_gene,"MANOVA result:", manova_pval) 
                #else: 
                    #print("skipped MANOVA because isoforms collinear")

            else: #if just one ILR column, us t.test / anova
                t = stats.ttest_ind([ilrs[0], ilrs[1], ilrs[2], ilrs[3], ilrs[4], ilrs[5], ilrs[6], ilrs[7], ilrs[8], ilrs[9], ilrs[10], ilrs[11]], [ilrs[12], ilrs[13], ilrs[14], ilrs[15], ilrs[16], ilrs[17], ilrs[18], ilrs[19], ilrs[20], ilrs[21], ilrs[22], ilrs[23]])
                print(new_gene, t) #print t-test output so we can tell which genes only have two isoforms
                print(new_gene, t.pvalue)
                #if t.pvalue < .05: #only print significant results
                #    print(new_gene, t)

