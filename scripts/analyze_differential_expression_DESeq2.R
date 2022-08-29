# Differential Expression analysis using DESeq2
# http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
BiocManager::install("DESeq2")
BiocManager::install("apeglm")
library(DESeq2)
library(tidyr)
library(tibble)
library(ggplot2)
library(patchwork)
library(apeglm)

#?DESeqDataSetFromHTSeqCount

sample_table_deseq2 <- read.table("analysis/DESeq2/sample_table.txt", header=T)
sample_table_deseq2$habitat <- as.factor(sample_table_deseq2$habitat)
counts <- read.table("analysis/DESeq2/htseq-count_out/htseq-count_results.2022-6-26.txt", row.names = 1) #with HAN412 gtf
#counts <- read.table("analysis/DESeq2/htseq-count_out/htseq-count_results.with_stringtie_gtf.2022-6-26.txt", row.names = 1) #with stringtie gtf
names(counts) <- paste(sample_table_deseq2$habitat,"_",sample_table_deseq2$sample_no, sep = "")

dds <- DESeqDataSetFromMatrix(countData = counts, #nrow(dds) = 46542
                              colData = sample_table_deseq2,
                              design = ~ habitat)

# pre-Filter out lowly expressed genes, where total read count is less than 24 across all samples combined
keep <- rowSums(counts(dds)) >= 24
dds_filtered <- dds[keep,]
## more stringent filter? only keep genes with minimum 240 reads across all samples, and with at least 10 samples having at least 10 reads. This is a very stringent filter, keeping approx only 50% of the data
#stringent_keep <- rowSums(counts(dds)) >= 240
#dds_filtered <- dds[stringent_keep,] #nrow(dds)=24370

#stringent_keep2 <- rowSums(counts(dds_filtered) >= 10) >= 10
#dds_filtered <- dds_filtered[stringent_keep2,] #nrow(dds) = 23803

# save list of all the expressed genes (>= 24 reads total) for use in GO Enrichment analysis as the 'background set' of genes
expressed_genes <- data.frame(Ha412_gene=rownames(DESeq2::counts(DE))) %>%
  mutate(Ha412_gene=gsub("mRNA","gene",Ha412_gene)) %>%
  mutate(Ha412_gene=gsub("rRNA","gene",Ha412_gene)) %>%
  mutate(Ha412_gene=gsub("ncRNA","gene",Ha412_gene)) %>%
  filter(grepl("gene", Ha412_gene))

write.table(expressed_genes, file = "analysis/DESeq2/expressed_genes.txt", quote = F, row.names = F, col.names = F, sep = "\t") #write list of all expressed genes. there are three rows at bottom of this list that need to be deleted, so total num of expressed genes is 32308

# get a list of all the lowly expressed genes that we filtered out, so we can filter these genes from our DEXSeq flattened gff file as well.
low_expressed_genes <- data.frame(counts(subset(dds,!(rownames(dds) %in% rownames(dds_filtered))))) %>%
  rownames_to_column("geneID") %>%
  mutate(geneID=gsub("mRNA:","gene_",geneID)) %>%
  dplyr::select(geneID)
write.table(low_expressed_genes, "low_expressed_genes.txt", quote = F, row.names = F, col.names = F)

# re-level the habitat factor to make non_dune samples the reference (i.e. the denominator in log Fold change calculation, log2(dune/non-dune))
dds_filtered$habitat <- relevel(dds_filtered$habitat, ref = "non-dune")

# Analyze differential expression
DE <- DESeq(dds_filtered)
save(DE, file="analysis/DESeq2/DESeq2_results.Robj") #save so we don't have to repeat this step.
load("analysis/DESeq2/DESeq2_results.Robj") #without stringtie
#load("analysis/DESeq2/DESeq2_results.stringtie.Robj") #with stringtie

# Extract normalized counts. Use these for eQTL analysis?
normalized_counts <- as.data.frame(DESeq2::counts(DE, normalized=TRUE))

# rlog transformation. for expression PCA. needs to be applied WITHOUT normalization bc it requires expression values to be integers
rlog_counts <- rlog(DESeq2::counts(DE, normalized=FALSE)) #32311 genes with <24 filter
# save to file for use in PCA
write.csv(rlog_counts, file = "analysis/pca/rlog-transformed_expression_counts.csv", quote = F, row.names = F)

temp <- data.frame(rlog_counts) %>% tibble::rownames_to_column("Ha412_gene")
write.table(temp, file="analysis/DESeq2/rlog-transformed_expression_counts.tsv", quote = F, row.names = F, sep = "\t")

#### DE results ####
# Provide contrast argument to set dune habitat in the numerator. Set LFC threshold of 1, which means a particular having twice as much expression in one ecotype vs the other. ALso, set significance threshold to .05 
#de_results <- data.frame(results(DE, contrast = c("habitat","dune","non-dune"), alpha = .05, cooksCutoff = T)) #%>% #no LFC threshold, no LFC shrinkage
#  rownames_to_column("id")
#de_results_p05 <- de_results %>%
#  arrange(desc(abs(log2FoldChange))) %>%
#  subset(padj<.05) #filter for just the significant DE genes
#de_results_lfc1 <- data.frame(results(DE, contrast = c("habitat","dune","non-dune"), alpha = .05, lfcThreshold #= 1, altHypothesis="greaterAbs", cooksCutoff = T)) %>% #threshold of LFC=1
#  rownames_to_column("id")
#de_results_lfc1_p05 <- data.frame(de_results_lfc1) %>%
#  arrange(desc(abs(log2FoldChange))) %>%
#  subset(padj<.05) #filter for just the significant DE genes (at LFC>=1)

# without logFold threshold.Shrink Log Fold Change for visualization and ranking
de_results_Shrink <- lfcShrink(DE, coef=c("habitat_dune_vs_non.dune"), type="apeglm")
summary(de_results_Shrink, alpha=.05)
de_results_p05_Shrink <- data.frame(de_results_Shrink) %>% # without stringtie
  rownames_to_column("Ha412_gene") %>% 
  mutate(Ha412_gene=gsub("mRNA","gene",Ha412_gene)) %>%
  mutate(Ha412_gene=gsub("ncRNA","gene",Ha412_gene)) %>%
  subset(padj<.05) %>% 
  filter(!grepl("Chr00", Ha412_gene), grepl("gene:", Ha412_gene)) %>% #exclude genes on non-chrom contigs and get rid of "non feature" and multi aln row.
  arrange(desc(abs(log2FoldChange))) %>%
  left_join(Ha412_Ath_mappings)
de_genes_noLFCthreshold <- subset(de_results_p05_Shrink, padj<.05) %>%
  dplyr::select(Ha412_gene, Ath_gene)
de_genes_dune_noLFCthreshold <- subset(de_results_p05_Shrink, log2FoldChange>0)
de_genes_non.dune_noLFCthreshold <- subset(de_results_p05_Shrink, log2FoldChange<0)

# with LFC=0.5 threshold
de_results_lfc0.5_Shrink <- lfcShrink(DE, coef=c("habitat_dune_vs_non.dune"), type="apeglm", lfcThreshold = .5)
summary(de_results_lfc0.5_Shrink) # 367 genes @ s<.005
de_results_lfc0.5_s005_Shrink <- data.frame(de_results_lfc0.5_Shrink) %>% # without stringtie
  rownames_to_column("Ha412_gene") %>% 
  mutate(Ha412_gene=gsub("mRNA","gene",Ha412_gene)) %>%
  mutate(Ha412_gene=gsub("ncRNA","gene",Ha412_gene)) %>%
  subset(svalue<.005) %>%
  filter(!grepl("Chr00", Ha412_gene)) %>%
  arrange(desc(abs(log2FoldChange))) %>%
  left_join(Ha412_Ath_mappings)
de_genes_dune_lfc0.5 <- subset(de_results_lfc0.5_s005_Shrink, log2FoldChange>0)
de_genes_non.dune_lfc0.5 <- subset(de_results_lfc0.5_s005_Shrink, log2FoldChange<0)

# with LFC=1 threshold. 
de_results_lfc1_Shrink <- lfcShrink(DE, coef=c("habitat_dune_vs_non.dune"), type="apeglm", lfcThreshold = 1)
summary(de_results_lfc1_Shrink) # 367 genes @ s<.005
de_results_lfc1_s005_Shrink <- data.frame(de_results_lfc1_Shrink) %>% # without stringtie
  rownames_to_column("Ha412_gene") %>% 
  mutate(Ha412_gene=gsub("mRNA","gene",Ha412_gene)) %>%
  mutate(Ha412_gene=gsub("ncRNA","gene",Ha412_gene)) %>%
  subset(svalue<.005) %>% 
  filter(!grepl("Chr00", Ha412_gene)) %>%
  arrange(desc(abs(log2FoldChange))) %>%
  left_join(Ha412_Ath_mappings)
de_genes_dune_lfc1 <- subset(de_results_lfc1_s005_Shrink, log2FoldChange>0)
de_genes_non.dune_lfc1 <- subset(de_results_lfc1_s005_Shrink, log2FoldChange<0)



# write list of genes to file
write.table(de_results_p05_Shrink, file = "analysis/DESeq2/DE_genes_noLFCthreshold_p05_Shrink.tsv", sep = '\t', quote = F, row.names = F)
write.table(de_results_p05_Shrink$Ha412_gene, file = "analysis/GO_analysis/study_DE_genes_noLFCthreshold.txt", sep = '\t', quote = F, row.names = F, col.names = F)
write.table(de_genes_dune_noLFCthreshold$Ha412_gene, file = "analysis/GO_analysis/study_DE_genes_dune_noLFCthreshold.txt", sep = '\t', quote = F, row.names = F, col.names = F)
write.table(de_genes_non.dune_noLFCthreshold$Ha412_gene, file = "analysis/GO_analysis/study_DE_genes_non-dune_noLFCthreshold.txt", sep = '\t', quote = F, row.names = F, col.names = F)

write.table(de_results_lfc0.5_s005_Shrink, file = "analysis/DESeq2/DE_genes_LFC0.5_sval005_Shrink.tsv", sep = '\t', quote = F, row.names = F)
write.table(de_genes_dune_lfc0.5$Ha412_gene, file = "analysis/GO_analysis/study_DE_genes_dune_LFC0.5.txt", sep = '\t', quote = F, row.names = F, col.names = F)
write.table(de_genes_non.dune_lfc0.5$Ha412_gene, file = "analysis/GO_analysis/study_DE_genes_non-dune_LFC0.5.txt", sep = '\t', quote = F, row.names = F, col.names = F)

write.table(de_results_lfc1_s005_Shrink, file = "analysis/DESeq2/DE_genes_LFC1_sval005_Shrink.tsv", sep = '\t', quote = F, row.names = F)
write.table(de_genes_dune$Ha412_gene, file = "analysis/GO_analysis/study_DE_genes_dune_LFC1.txt", sep = '\t', quote = F, row.names = F, col.names = F)
write.table(de_genes_non.dune$Ha412_gene, file = "analysis/GO_analysis/study_DE_genes_non-dune_LFC1.txt", sep = '\t', quote = F, row.names = F, col.names = F)

#write.table(dplyr::select(de_results_p05, stringtie_gene), file = "analysis/DESeq2/DE_genes_padj.05.tsv", sep = '\t', quote = F, row.names = F)
#write.table(dplyr::select(de_results_lfc1_p05, stringtie_gene), file = "analysis/DESeq2/DE_genes_LFC1_padj.05.tsv", sep = '\t', quote = F, row.names = F)


#### Map differential expression across the genome ####
# Read-in our annotations, for the gene positions. 
gff <- data.frame(rtracklayer::import("data/ref_genome_Ha412HO/HAN412_Eugene_curated_v1_1.gff3")) %>%
  rename(seqnames="chrom", Parent="Ha412_gene")
# Filter for only the mRNA, ncRNA, and rRNA (has same positions as 'gene'; there are 8 ncRNAs in or DE gene list, notably)
transcripts <- filter(gff, type %in% c("mRNA", "ncRNA", "rRNA"))

## use stringtie gtf instead
#gtf <- data.frame(rtracklayer::import("data/stringtie/merged_transcripts.stringtie.gtf")) %>%
#  rename(seqnames="chrom") 
#transcripts <- subset(gtf, type=="transcript")

# get the cumulative position of each gene and then the max gene end position on each chrom
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

# expression manhattan plot but with no LFC threshold on DE test
expr_manhattan_df_noLFCthreshold <- left_join(data.frame(de_results_Shrink) %>% rownames_to_column("ID"), dplyr::select(transcripts, ID, chrom, start, end, strand)) %>%
  inner_join(cum_lengths) %>% #inner join gets rid of non-chromosomal contigs, So we lose a few DE genes here.
  mutate(startcum = start + tot, sig=if_else(padj<.05, "yes", "no"))

sig_data_expr <- expr_manhattan_df_noLFCthreshold %>% 
  subset(padj < .05)
notsig_data_expr <- expr_manhattan_df_noLFCthreshold %>% 
  subset(padj >= .05) %>%
  group_by(chrom) %>% 
  sample_frac(0.75)
notsig_remove_expr <- subset(notsig_data_expr,abs(log2FoldChange) < .05) %>%
  sample_frac(0.95)
notsig_keep_expr <- !(notsig_data_expr$ID %in% notsig_remove_expr$ID)
notsig_data_expr <- notsig_data_expr[notsig_keep_expr,]

expr_manhattan_df_noLFCthreshold_reduced <- bind_rows(sig_data_expr, notsig_data_expr)

# get middle position of each chromosome for x axis label
axisdf_expr <- expr_manhattan_df_noLFCthreshold_reduced %>%
  group_by(chrom) %>%
  summarize(center=(max(startcum) + min(startcum)) / 2) %>%
  mutate(chrom=c("1","2","3","4","5",
                 "6","7","8","9","10",
                 "11","12","13","14","15",
                 "16","17"))

expr_manhattan_plot <- ggplot(expr_manhattan_df_noLFCthreshold_reduced, aes(x=startcum, y=log2FoldChange, color=as.factor(chrom), alpha=as.factor(sig)=="yes")) +
  geom_point(size=3) +
  scale_x_continuous(label = axisdf_expr$chrom, breaks = axisdf_expr$center, expand=c(0.0175,0.0175)) +
  scale_color_manual(values = rep(c("cornflowerblue", "grey50"), unique(length(axisdf_expr$chrom)))) +
  scale_alpha_manual(values = c(.1,1)) +
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(size = 36),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  #geom_hline(yintercept = -1, linetype = "dashed") +
  #geom_hline(yintercept = 1, linetype = "dashed") +
  geom_segment(data=inversion_regions, aes(x=startcum, xend=endcum, y=0, yend=0),
               size=6, color="red", alpha=0.8) +
  labs(x="", title = "Differential expression")

ggsave("figures/expr_manhattan.png",plot = expr_manhattan_plot, device = "png", width = 8, height = 6, units = "in", dpi = 300)

#expr_manhattan_df <- left_join(data.frame(de_results_lfc1_Shrink) %>% rownames_to_column("ID"), dplyr::select(transcripts, ID, #chrom, start, end, strand)) %>%
#  inner_join(cum_lengths) %>% #inner join gets rid of non-chromosomal contigs. There were 12 sig DE genes on these unmapped contigs#. 
#  #rowwise() %>%
#  # get the cumulative start position of each gene
#  # and yes/no if it passes significance threshold, in this case svalue<.005
#  mutate(startcum = start + tot, sig=if_else(svalue<.005, "yes", "no"))
#
## split apart significant DE genes and nonsignifcant DE genes so we can downsample the nonsig data #for prettier plot
#sig_data_expr <- expr_manhattan_df %>% 
#  subset(svalue < .005)
#notsig_data_expr <- expr_manhattan_df %>% 
#  subset(svalue >= .005) %>%
#  group_by(chrom) %>% 
#  sample_frac(0.75)
#notsig_remove_expr <- subset(notsig_data,abs(log2FoldChange) < .05) %>%
#  sample_frac(0.95)
#notsig_keep_expr <- !(notsig_data$ID %in% notsig_remove$ID)
#notsig_data_expr <- notsig_data[notsig_keep,]
#
#expr_manhattan_df_reduced <- bind_rows(sig_data_expr, notsig_data_expr)
#
## get middle position of each chromosome for x axis label
#axisdf_expr <- expr_manhattan_df_reduced %>%
#  group_by(chrom) %>%
#  summarize(center=(max(startcum) + min(startcum)) / 2) %>%
#  mutate(chrom=c("1","2","3","4","5",
#                 "6","7","8","9","10",
#                 "11","12","13","14","15",
#                 "16","17"))
#
#expr_manhattan_plot <- ggplot(expr_manhattan_df_reduced, aes(x=startcum, y=log2FoldChange, color=as.factor(chrom), alpha=as.factor#(sig)=="yes")) +
#  geom_point(size=3) +
#  scale_x_continuous(label = axisdf_expr$chrom, breaks = axisdf_expr$center, expand=c(0.0175,0.0175)) +
#  scale_color_manual(values = rep(c("cornflowerblue", "grey50"), unique(length(axisdf_expr$chrom)))) +
#  scale_alpha_manual(values = c(.1,1)) +
#  theme_classic() +
#  theme(legend.position = "none",
#        text = element_text(size = 36),
#        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
#  geom_hline(yintercept = -1, linetype = "dashed") +
#  geom_hline(yintercept = 1, linetype = "dashed") +
#  geom_segment(data=inversion_regions, aes(x=startcum, xend=endcum, y=0, yend=0),
#               size=6, color="red", alpha=0.8) +
#  labs(x="Chromosome", title = "Differential expression")
#ggsave("figures/expr_manhattan.png",plot = expr_manhattan_plot, device = "png", width = 8, height = 6, units = "in", dpi = 300)

#### how many DE genes on each chromosome? ####
chrom_lengths <- read.table("data/ref_genome_Ha412HO/chrom_sizes_Ha412HOv2.0.txt", col.names = c("chrom", "length"))
de_genes_per_chrom <- de_results %>% left_join(dplyr::select(transcripts, ID, chrom)) %>%
  subset(padj<.05) %>%
  group_by(chrom) %>%
  filter(!grepl("Chr00", chrom)) %>%
  summarize(num_de_genes=n()) %>%
  right_join(chrom_lengths) %>%
  mutate(num_de_genes_norm = (num_de_genes/length)*1000000)
de_genes_per_chrom$chrom <- c("01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17")

DE_genes_per_chrom_plot <- ggplot(data=de_genes_per_chrom, aes(x=chrom, y=num_de_genes_norm, fill=chrom)) +
  geom_bar(stat = "identity", alpha=0.9) +
  theme_classic() +
  scale_fill_manual(values = rep(c("cornflowerblue", "grey50"), unique(length(DE_genes_per_chrom_plot)))) +
  labs(x="Chromosome", y="DE genes per 1Mb") +
  theme(text = element_text(size=36),
        legend.position = "none",
        legend.title = element_blank())

#FIG_2 <- plot_spacer() / Fst_manhattan_plot | (expr_manhattan_plot / DE_genes_per_chrom_plot)
#ggsave("figures/FIG_2.jpg", plot=FIG_2, width=56, height = 28, dpi = 600, units = "cm", bg = "white")


#### MA plot ####
MA_plot <- ggplot(data.frame(de_results_lfc1_Shrink), aes(x=baseMean, y=log2FoldChange, color=svalue<.005)) +
  geom_point(alpha=0.5, size=3, show.legend = F) +
  scale_color_manual(values = c('TRUE'='cornflowerblue', 'FALSE'='grey50')) +
  scale_x_log10() +
  theme_classic()

DESeq2::plotMA(de_results_Shrink)
DESeq2::plotMA(DE)
