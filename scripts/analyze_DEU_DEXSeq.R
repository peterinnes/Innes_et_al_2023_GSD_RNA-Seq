BiocManager::install("DEXSeq")
BiocManager::install("BiocParallel")
library(DEXSeq)
library(BiocParallel)
library(gtools)
library(dplyr)

# see vignette: https://rdrr.io/bioc/DEXSeq/f/vignettes/DEXSeq.Rmd

# read-in the data
count_files <- list.files("analysis/DEXSeq/dexseq-count_out/stringtie", pattern = ".dexseq_counts.filtered.txt", full.names = T) #with stringtie, counted after pre-filter of low-expressed genes from dexseq-flattened gtf file
#count_files <- list.files("analysis/DEXSeq/dexseq-count_out", pattern = ".dexseq_counts.noaggreg.filtered.txt", full.names = T) #without stringtie, counted after pre-filter of low-expressed genes from dexseq-flattened gtf file
count_files <- mixedsort(count_files, decreasing = T) #reorder with gtools::mixedsort() to match the order in sample_table_dexseq. not sure if this matters.

# use the pre-filtered annotations (and thus filtered counts, above) per Soneson et al 2016. We exclude (pre-filter) the same set of lowly expressed genes as identified in the pre-filtering step of DESeq2 pipeline: genes with fewer than 24 reads total.
annotations <- list.files("data/ref_genome_Ha412HO", pattern = "DEXSeq.filtered.gff", full.names = T) #without stringtie
#annotations <- list.files("data/stringtie", pattern = "DEXSeq.filtered.gff", full.names = T) #with stringtie
sample_table_dexseq <- read.table("sample_list.txt", header=T)
rownames(sample_table_dexseq) <- sample_table_dexseq$sample
sample_table_dexseq <- sample_table_dexseq[1]
sample_table_dexseq$habitat <- as.factor(sample_table_dexseq$habitat)

# construct a DEXSeqDataSet object
dxd <- DEXSeqDataSetFromHTSeq(count_files, sampleData = sample_table_dexseq,
                              design = ~ sample + exon + habitat:exon,
                              flattenedfile = annotations)

# normalise and estimate dispersion
BPPARAM <- MulticoreParam(8)
dxd <- estimateSizeFactors(dxd)
dxd <- estimateDispersions(dxd, BPPARAM = BPPARAM)

# test for DEU
dxd <- testForDEU(dxd, BPPARAM = BPPARAM)
dxd <- estimateExonFoldChanges( dxd, fitExpToVar="habitat", denominator = "non-dune", BPPARAM = BPPARAM) #estimate relative exon usage fold changes

# results
dxr <- DEXSeqResults( dxd )
save(dxr, file="analysis/DEXSeq/DEXSeq_results.filtered.Robj") #save so we don't have to repeat this step
load("analysis/DEXSeq/DEXSeq_results.filtered.Robj")

dxr_df <- data.frame(dxr) #convert to dataframe
#length(unique(dxr_df$groupID)) #how many genes total?

# how many DEU exons? 1425 @ alpha=.05, without pre-filter. 1774 with filter. 3099 wtih stringtie gtf and prefilter.
table(dxr$padj < .05)

# How many DEU genes? Need to control for multiple testing at the gene level. 
num_deu_genes <- sum( perGeneQValue(dxr) < 0.05)
num_deu_genes #956 with pre-filtered, noaggreg dataset; 452 with stringtie gtf-aligned, pre-filtered, noaggreg data set.

deu_genes_noLFCthreshold <- data.frame(perGeneQValue(dxr)) %>%
  rownames_to_column("groupID") %>%
  subset(perGeneQValue.dxr. < .05) %>%
  mutate(groupID = gsub("gene_", "gene:", groupID)) %>%
  dplyr::select(groupID) %>%
  rename(groupID="Ha412_gene") %>%
  left_join(Ha412_Ath_mappings) #match with Ath gene names

# genes with at least one exon passing threshold of LFC=1, to be used in a subsequent filtering step.
deu_list_lfc1 <- dxr_df %>%
  subset(abs(log2fold_dune_non.dune)>=1) %>%
  mutate(groupID = gsub("gene_", "", groupID)) %>%
  dplyr::select(groupID) %>%
  unique()

# get list of the DEU genes at fdr of .05 and at least one exon @ LFC=1. 448 without stringtie, with noaggreg, prefiltered. 354 genes, with the stringtie gtf-aligned, no aggreg, and pre-filtering. 
deu_genes_lfc1 <- data.frame(perGeneQValue(dxr)) %>%
  rownames_to_column("groupID") %>%
  subset(perGeneQValue.dxr. < .05) %>%
  mutate(groupID = gsub("gene_", "", groupID)) %>%
  dplyr::select(groupID) %>%
  subset(groupID %in% deu_list_lfc1$groupID) %>% #filter for LFC=1 threshold
  mutate(groupID=gsub("Ha412","gene:Ha412", groupID)) %>% #add 'gene' so we can join with Ha412 gene names and Ath gene names
  rename(groupID="Ha412_gene") %>%
  left_join(Ha412_Ath_mappings) #match with Ath gene names
rownames(deu_genes_lfc1) <- NULL

write.table(deu_genes_noLFCthreshold, file="analysis/DEXSeq/DEU_genes_padj.05.tsv", row.names = F, quote = F, sep = '\t')
write.table(deu_genes_noLFCthreshold$Ha412_gene, file="analysis/GO_analysis/study_DEU_genes_noLFCthreshold_padj.05.txt", row.names = F, quote = F, sep = '\t', col.names = F)
write.table(deu_genes_lfc1, file="analysis/DEXSeq/DEU_genes_LFC1_padj.05.tsv", row.names = F, quote = F, sep = '\t')
# plot DEU for Ha412HOChr09g0373721, the only gene found to be DS in all three splicing analyses (DEXSeq, rMATS, and parents_diff_v2.py)
plotDEXSeq(object = dxr, geneID = "gene_Ha412HOChr05g0235401", splicing = T, expression = F, fitExpToVar="habitat",displayTranscripts = F, legend = T, color = c("gold2", "forestgreen"),)


#### make dataframe of rlog-transforms exon counts for use in PCA ####
rlog_exon_counts <- DESeq2::rlog(DEXSeq::counts(dxr, normalized=FALSE))
write.csv(rlog_exon_counts, file = 
            "analysis/pca/dexseq_exon_counts.csv", quote = F, row.names = F)


#### plot DEU across genome ####
chrom_lengths <- read.table("Ha412H0_chrom_lengths.txt", col.names = c("chrom", "length")) %>%
  mutate(length=as.numeric(length))
cum_lengths <- chrom_lengths %>%
  mutate(tot=cumsum(length)-length) %>%
  dplyr::select(-length)

# get inversion regions
inversion_regions <- read.table("analysis/inversions/inversion_regions.txt") %>%
  dplyr::select(V2,V3,V4,V9) %>%
  rename(V2="chrom", V3="start", V4="end", V9="name") %>%
  left_join(cum_lengths) %>%
  mutate(startcum=start+tot, endcum=end+tot)

dxr_df$transcripts <- unlist(dxr_df$transcripts)
#absmax <- function(x) { x[which.max( abs(x) )]} #function to get maximum signed LFC

# subset our DEU dataframe for the exon with max(abs(log2FoldChange))
dxr_maxLFC_df <- data.frame(dplyr::select(dxr_df, genomicData.seqnames, genomicData.start, genomicData.end, transcripts, groupID,featureID, log2fold_dune_non.dune) %>%
  rename(genomicData.seqnames = "chrom", genomicData.start="start", genomicData.end="end", transcripts="ID") %>%
  #na.omit() %>%
  group_by(ID) %>%
  filter(n_distinct(featureID)>1) %>%
  arrange(desc(abs(log2fold_dune_non.dune))) %>%
  dplyr::slice(1) %>%
  ungroup())

deu_manhattan_df <- inner_join(dxr_maxLFC_df, cum_lengths) %>% #inner join gets rid of non-chromosomal contigs.
  left_join(data.frame(perGeneQValue(dxr)) %>% rownames_to_column("groupID")) %>%
  mutate(startcum = start + tot, sig=if_else(perGeneQValue.dxr.<.05, "yes", "no"))

axisdf_deu <- deu_manhattan_df %>%
  group_by(chrom) %>%
  summarize(center=(max(startcum) + min(startcum)) / 2) %>%
  mutate(chrom=c("1","2","3","4","5",
                 "6","7","8","9","10",
                 "11","12","13","14","15",
                 "16","17"))

sig_deu_data <- deu_manhattan_df %>% 
  subset(sig=="yes")
notsig_deu_data <- deu_manhattan_df %>% 
  subset(sig=="no") %>%
  group_by(chrom) %>% 
  sample_frac(0.5)

deu_manhattan_df_reduced <- bind_rows(sig_deu_data, notsig_deu_data)
deu_manhattan_plot <- ggplot(deu_manhattan_df_reduced, aes(x=startcum, y=log2fold_dune_non.dune, color=as.factor(chrom), alpha=as.factor(sig)=="yes")) +
  geom_point(size=3) +
  scale_x_continuous(label = axisdf_deu$chrom, breaks = axisdf_deu$center, expand=c(0.0175,0.0175)) +
  scale_y_continuous(limits=c(-21,20)) +
  scale_color_manual(values = rep(c("cornflowerblue", "grey50"), unique(length(axisdf_deu$chrom)))) +
  scale_alpha_manual(values = c(.1,1)) +
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(size = 36),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  #geom_hline(yintercept = -1, linetype = "dashed") +
  #geom_hline(yintercept = 1, linetype = "dashed") +
  #geom_segment(data=inversion_regions_expr, aes(x=startcum, xend=endcum, y=0, yend=0),
               #size=6, color="red", alpha=0.8) +
  labs(x="Chromosome", y="rel. Exon Usage LFC", title = "Differential splicing\n(exon-based)") +
  geom_segment(data=inversion_regions, aes(x=startcum, xend=endcum, y=0, yend=0),
               size=6, color="red", alpha=0.8)

ggsave("figures/deu_manhattan.png",plot = deu_manhattan_plot, device = "png", width = 8, height = 6, units = "in", dpi = 300)
