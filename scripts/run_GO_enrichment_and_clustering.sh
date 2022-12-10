# this script needs to be run with 'bash -i' so that conda can be activated and deactivated properly
cd ~/gsd_RNA-seq/analysis/GO_analysis

# run Gene Ontology enrichment analyses for our various gene sets of interest
conda activate goatools

find_enrichment.py --obo go-basic.obo study_DE_genes_dune_noLFCthreshold.txt expressed_genes.txt Ha412_GO_GOSLIM.id2gos.txt \
	--pval=.05 --method fdr_bh --pval_field=fdr_bh --no_propagate_counts \
	--outfile=results_DE_dune_noLFCthreshold_GO_enrichment.txt | tee stderr_GO_enrichment.DE_dune_noLFCthreshold.txt
awk '$3!="p"' results_DE_dune_noLFCthreshold_GO_enrichment.txt > tmp
mv tmp results_DE_dune_noLFCthreshold_GO_enrichment.txt

find_enrichment.py --obo go-basic.obo study_DE_genes_non-dune_noLFCthreshold.txt expressed_genes.txt Ha412_GO_GOSLIM.id2gos.txt \
	--pval=.05 --method fdr_bh --pval_field=fdr_bh --no_propagate_counts \
	--outfile=results_DE_non-dune_noLFCthreshold_GO_enrichment.txt | tee stderr_GO_enrichment.DE_non-dune_noLFCthreshold.txt
awk '$3!="p"' results_DE_non-dune_noLFCthreshold_GO_enrichment.txt > tmp
mv tmp results_DE_non-dune_noLFCthreshold_GO_enrichment.txt

find_enrichment.py --obo go-basic.obo study_DS_rMATS_genes.txt expressed_multi_exonic_genes.txt Ha412_GO_GOSLIM.id2gos.txt \
	--pval=.05 --method fdr_bh --pval_field=fdr_bh --no_propagate_counts \
	--outfile=results_DS_rMATS_GO_enrichment.txt | tee stderr_GO_enrichment.DS_rMATS.txt
awk '$3!="p"' results_DS_rMATS_GO_enrichment.txt > tmp
mv tmp results_DS_rMATS_GO_enrichment.txt

find_enrichment.py --obo go-basic.obo study_intersect_DE–DS_genes.txt expressed_multi_exonic_genes.txt Ha412_GO_GOSLIM.id2gos.txt \
	--pval=.05 --method fdr_bh --pval_field=fdr_bh --no_propagate_counts \
	--outfile=results_intersect_DE–DS_noLFCthreshold_GO_enrichment.txt | tee stderr_GO_enrichment.intersect_DE–DS_noLFCthreshold.txt
awk '$3!="p"' results_intersect_DE–DS_noLFCthreshold_GO_enrichment.txt > tmp
mv tmp results_intersect_DE–DS_noLFCthreshold_GO_enrichment.txt

conda deactivate

# cluster GO terms based on similarity and overlap in genes, using GOMCL
mkdir GOMCL_BP GOMCL_CC GOMCL_MF

# cluster all results except for the DE–DS intersect gene set bc it's so short.
# have to run the clustering seprately for each GO category (biological process, cellular component, molecular function)
GOMCL.py -got GOATOOLS -gotype BP -nw go-basic.obo results_DE_dune_noLFCthreshold_GO_enrichment.txt
GOMCL.py -got GOATOOLS -gotype BP -nw go-basic.obo results_DE_non-dune_noLFCthreshold_GO_enrichment.txt
GOMCL.py -got GOATOOLS -gotype BP -nw go-basic.obo results_DS_rMATS_GO_enrichment.txt
mv results_*GOsize* GOMCL_BP/

GOMCL.py -got GOATOOLS -gotype CC -nw go-basic.obo results_DE_dune_noLFCthreshold_GO_enrichment.txt
GOMCL.py -got GOATOOLS -gotype CC -nw go-basic.obo results_DE_non-dune_noLFCthreshold_GO_enrichment.txt
GOMCL.py -got GOATOOLS -gotype CC -nw go-basic.obo results_DS_rMATS_GO_enrichment.txt
mv results_*GOsize* GOMCL_CC/

GOMCL.py -got GOATOOLS -gotype MF -nw go-basic.obo results_DE_dune_noLFCthreshold_GO_enrichment.txt
GOMCL.py -got GOATOOLS -gotype MF -nw go-basic.obo results_DE_non-dune_noLFCthreshold_GO_enrichment.txt
GOMCL.py -got GOATOOLS -gotype MF -nw go-basic.obo results_DS_rMATS_GO_enrichment.txt
mv results_*GOsize* GOMCL_MF/

