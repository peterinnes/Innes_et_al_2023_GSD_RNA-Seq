BiocManager::install("DEXSeq")
BiocManager::install("BiocParallel")
library(DEXSeq)
library(BiocParallel)
library(gtools)
library(dplyr)

# see vignette: https://rdrr.io/bioc/DEXSeq/f/vignettes/DEXSeq.Rmd

# read-in the data
count_files <- list.files("analysis/DEXSeq/dexseq-count_out", pattern = ".dexseq_counts.txt", full.names = T)
count_files <- mixedsort(count_files, decreasing = T) #reorder wtih gtools::mixedsort() to match the sample_table. not sure if this matters.

annotations <- list.files("data/ref_genome_Ha412HO", pattern = "DEXSeq.gff", full.names = T)
sample_table <- read.table("sample_list.txt", header=T)
rownames(sample_table) <- sample_table$sample
sample_table <- sample_table[1]
sample_table$habitat <- as.factor(sample_table$habitat)

# construct a DEXSeqDataSet object
dxd <- DEXSeqDataSetFromHTSeq(count_files, sampleData = sample_table,
                              design = ~ sample + exon + habitat:exon,
                              flattenedfile = annotations)

# normalise and estimate dispersion
BPPARAM <- MulticoreParam(4)
dxd <- estimateSizeFactors(dxd)
dxd <- estimateDispersions(dxd, BPPARAM = BPPARAM)

# test for DEU
dxd <- testForDEU(dxd, BPPARAM = BPPARAM)
dxd <- estimateExonFoldChanges( dxd, fitExpToVar="habitat", denominator = "non-dune", BPPARAM = BPPARAM) #estimate relative exon usage fold changes

# results
dxr <- DEXSeqResults( dxd )
save(dxr, file="analysis/DEXSeq/DEXSeq_results.Robj") #save so we don't have to repeat this step
load("analysis/DEXSeq/DEXSeq_results.Robj")

dxr_df <- data.frame(dxr) #convert to dataframe
#length(unique(dxr_df$groupID)) #how many genes total?

# genes with at least one exon passing threshold of LFC=1 (disregarding signficance here)
lfc1_deu_genes <- dxr_df %>%
  subset(abs(log2fold_dune_non.dune)>=1) %>%
  mutate(groupID = gsub("gene_", "", groupID)) %>%
  dplyr::select(groupID) %>%
  unique()

# how many DEU exons? 1425 @ alpha=.05
table(dxr$padj < .05)

# How many DEU genes? Need to control for multiple testing at the gene level. 
num_deu_genes <- sum( perGeneQValue(dxr) < 0.05)
num_deu_genes #789 genes

# get list of the DEU genes at fdr of .05 and at least one exon @ LFC=1
deu_genes <- data.frame(perGeneQValue(dxr)) %>%
  rownames_to_column("groupID") %>%
  subset(perGeneQValue.dxr. < .05) %>%
  mutate(groupID = gsub("gene_", "", groupID)) %>%
  dplyr::select(groupID) %>%
  subset(!groupID %in% low_expressed_genes$groupID) %>% #filter out lowly expressed genes
  #subset(groupID %in% lfc1_deu_genes$groupID) %>% #filter for LFC=1 threshold
  mutate(groupID=gsub("Ha412","mRNA:Ha412", groupID)) #add 'mRNA:' label so that it's in usable format for the GO analysis
rownames(deu_genes) <- NULL

write.table(deu_genes, file="analysis/DEXSeq/DEU_genes_padj.05.tsv", row.names = F, quote = F, sep = '\t')

deu_genes <- deu_genes %>%
  rename(groupID="Ha412_gene") %>%
  left_join(Ha412_Ath_mappings)
