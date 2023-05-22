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
library(rtracklayer)

#?DESeqDataSetFromHTSeqCount

sample_table_deseq2 <- read.table("analysis/DESeq2/sample_table.txt", header=T)
sample_table_deseq2$habitat <- as.factor(sample_table_deseq2$habitat)
counts <- read.table("analysis/DESeq2/htseq-count_out/htseq-count_results.2022-6-26.txt",
                     row.names = 1)


names(counts) <- paste(sample_table_deseq2$habitat,"_",
                       sample_table_deseq2$sample_no, sep = "")

dds <- DESeqDataSetFromMatrix(countData = counts, #nrow(dds) = 46542
                              colData = sample_table_deseq2,
                              design = ~ habitat)

# pre-Filter out lowly expressed genes,
# where read count is less than 24 across all samples combined
keep <- rowSums(counts(dds)) >= 24
dds_filtered <- dds[keep,]

## more stringent filter? only keep genes with minimum 240 reads
## across all samples, and with at least 10 samples
## having at least 10 reads. This is a very stringent filter,
## keeping approx only 50% of the data
#stringent_keep <- rowSums(counts(dds)) >= 240
#dds_filtered <- dds[stringent_keep,] #nrow(dds)=24370
#stringent_keep2 <- rowSums(counts(dds_filtered) >= 10) >= 10
#dds_filtered <- dds_filtered[stringent_keep2,] #nrow(dds) = 23803

# save list of all the expressed genes (>= 24 reads cumendal)
# for use in GO Enrichment analysis as the 'background set' of genes
expressed_genes <- data.frame(Ha412_gene=rownames(counts(dds_filtered))) %>%
  mutate(Ha412_gene=gsub(".*RNA","gene",Ha412_gene)) %>%
  filter(grepl("gene", Ha412_gene))

# there are three rows at bottom of this list that need to be deleted,
# so num of expressed genes is 32308
write.table(expressed_genes, file = "analysis/DESeq2/expressed_genes.txt",
            quote = F, row.names = F, col.names = F, sep = "\t") 

# get a list of all the lowly expressed genes that we filtered out,
# so we can filter these genes from our DEXSeq flattened gff file as well.
low_expressed_genes <- data.frame(counts(
  subset(dds,!(rownames(dds) %in% rownames(dds_filtered))))) %>%
  rownames_to_column("geneID") %>%
  mutate(geneID=gsub("mRNA:","gene_",geneID)) %>%
  dplyr::select(geneID)
write.table(low_expressed_genes, "low_expressed_genes.txt", quote = F,
            row.names = F, col.names = F)


# get list of all multi-exonic genes in the genome
gff <- data.frame(rtracklayer::
                    import("data/ref_genome_Ha412HO/HAN412_Eugene_curated_v1_1.gff3")) %>%
  rename(chrom=seqnames, Ha412_gene=Parent)

multi_exon_genes <- filter(gff, type=="exon") %>%
  group_by(Ha412_gene) %>%
  summarise(num_exons=n()) %>%
  mutate(Ha412_gene=unlist(Ha412_gene),
         Ha412_gene=gsub(".*RNA", "gene", Ha412_gene)) %>%
  filter(num_exons > 1)

expressed_multi_exonic_genes <- inner_join(multi_exon_genes, expressed_genes)
write.table(expressed_multi_exonic_genes[1], "analysis/GO_analysis/expressed_multi_exonic_genes.txt",
            row.names = F, quote = F, col.names = F)

# re-level the habitat factor to make non_dune samples the reference
# (i.e. the denominator in log Fold change calculation, log2(dune/non-dune))
dds_filtered$habitat <- relevel(dds_filtered$habitat, ref = "non-dune")

# Analyze differential expression
DE <- DESeq(dds_filtered)
save(DE, file="analysis/DESeq2/DESeq2_results.Robj")
load("analysis/DESeq2/DESeq2_results.Robj") #without stringtie
#load("analysis/DESeq2/DESeq2_results.stringtie.Robj") #with stringtie

# normalized counts will be used later for eQTL analysis?
normalized_counts <- as.data.frame(DESeq2::counts(DE, normalized=TRUE))
save(normalized_counts, file = "data2/Rdata/DESeq2_normalized_counts.Rdata")

# regularized-log transformation, for expression PCA. Needs to be applied
# WITHOUT normalization bc it requires expression values to be integers
rlog_counts <- rlog(DESeq2::counts(DE, normalized=FALSE)) #32308 genes w/ >=24 reads
write.csv(rlog_counts, file = "analysis/pca/
          rlog-transformed_expression_counts.csv", quote = F, row.names = F)

temp <- data.frame(rlog_counts) %>% tibble::rownames_to_column("Ha412_gene")
write.table(temp, file="analysis/DESeq2/rlog-transformed_expression_counts.tsv",
            quote = F, row.names = F, sep = "\t")

#### DE results ####
# DE results without logFold threshold. Shrink Log Fold Change for visualization and ranking
de_results_Shrink <- lfcShrink(DE, coef=c("habitat_dune_vs_non.dune"), type="apeglm")
summary(de_results_Shrink, alpha=.05)

## Provide contrast argument to set dune habitat in the numerator. Set LFC threshold of 1, which means a particular gene having twice as much expression in one ecotype vs the other. ALso, set significance threshold to .05 
#de_results_lfc1 <- data.frame(results(DE, contrast = c("habitat","dune","non-dune"), alpha = .05, lfcThreshold = 1, altHypothesis="greaterAbs", cooksCutoff = T)) %>% rownames_to_column("id")

de_results_Shrink_df <- data.frame(de_results_Shrink) %>% #table for Sup Mat
  rownames_to_column("Ha412_gene") %>%
  mutate(Ha412_gene=gsub(".*RNA","gene",Ha412_gene)) %>%
  filter(grepl("gene:", Ha412_gene)) %>%
  arrange(padj) %>%
  left_join(Ha412_Ath_mappings) %>%
  rename(Ha412_gene="Ha412HOv2_gene",
         Ath_gene="Araport11_gene")
  
de_results_p05_Shrink_df <- data.frame(de_results_Shrink) %>% # without stringtie
  rownames_to_column("Ha412_gene") %>% 
  mutate(Ha412_gene=gsub(".*RNA","gene",Ha412_gene)) %>% #replace 'mRNA' and 'ncRNA' w/ 'gene'
  subset(padj<.05) %>% 
  #filter(!grepl("Chr00", Ha412_gene)) %>% #exclude genes on non-chrom contigs
  filter(grepl("gene:", Ha412_gene)) %>% #get rid of 'non feature' and 'multi aln'
  #arrange(desc(abs(log2FoldChange))) %>%
  arrange(padj) %>%
  left_join(Ha412_Ath_mappings)

de_genes_dune_noLFCthreshold <- subset(de_results_p05_Shrink_df, log2FoldChange>0)
de_genes_non.dune_noLFCthreshold <- subset(de_results_p05_Shrink_df, log2FoldChange<0)

## with LFC=0.5 threshold
#de_results_lfc0.5_Shrink <- lfcShrink(DE, coef=c("habitat_dune_vs_non.dune"),
#                                      type="apeglm", lfcThreshold = .5)
#summary(de_results_lfc0.5_Shrink) # 367 genes @ s<.005
#de_results_lfc0.5_s005_Shrink <- data.frame(de_results_lfc0.5_Shrink) %>%
#  rownames_to_column("Ha412_gene") %>% 
#  mutate(Ha412_gene=gsub("mRNA","gene",Ha412_gene)) %>%
#  mutate(Ha412_gene=gsub("ncRNA","gene",Ha412_gene)) %>%
#  subset(svalue<.005) %>%
#  filter(!grepl("Chr00", Ha412_gene)) %>%
#  arrange(desc(abs(log2FoldChange))) %>%
#  left_join(Ha412_Ath_mappings)
#de_genes_dune_lfc0.5 <- subset(de_results_lfc0.5_s005_Shrink, log2FoldChange>0)
#de_genes_non.dune_lfc0.5 <- subset(de_results_lfc0.5_s005_Shrink, log2FoldChange<0)

## with LFC=1 threshold. 
#de_results_lfc1_Shrink <- lfcShrink(DE, coef=c("habitat_dune_vs_non.dune"),
#                                    type="apeglm", lfcThreshold = 1)
#summary(de_results_lfc1_Shrink) # 367 genes @ s<.005
#de_results_lfc1_s005_Shrink <- data.frame(de_results_lfc1_Shrink) %>% # without stringtie
#  rownames_to_column("Ha412_gene") %>% 
#  mutate(Ha412_gene=gsub("mRNA","gene",Ha412_gene)) %>%
#  mutate(Ha412_gene=gsub("ncRNA","gene",Ha412_gene)) %>%
#  subset(svalue<.005) %>% 
#  filter(!grepl("Chr00", Ha412_gene)) %>%
#  arrange(desc(abs(log2FoldChange))) %>%
#  left_join(Ha412_Ath_mappings)
#de_genes_dune_lfc1 <- subset(de_results_lfc1_s005_Shrink, log2FoldChange>0)
#de_genes_non.dune_lfc1 <- subset(de_results_lfc1_s005_Shrink, log2FoldChange<0)

# save DE results as Rdata and write tables/gene lists to file
# e.g. for GO analysis
save(de_results_Shrink, file = "data2/Rdata/DESeq2_results_Shrink.Rdata")
write.table(de_results_Shrink_df,
            file="analysis/DESeq2/DESeq2_results_shrinkage.tsv", sep = '\t',
            quote = F, row.names = F) # Extended Data S1
write.table(de_results_p05_Shrink_df,
            file = "analysis/DESeq2/DE_genes_noLFCthreshold_p05_Shrink.tsv",
            sep = '\t', quote = F, row.names = F)
write.table(de_results_p05_Shrink_df$Ha412_gene,
            file = "analysis/GO_analysis/study_DE_genes_noLFCthreshold.txt",
            sep = '\t', quote = F, row.names = F, col.names = F)
write.table(de_genes_dune_noLFCthreshold$Ha412_gene,
            file = "analysis/GO_analysis/study_DE_genes_dune_noLFCthreshold.txt",
            sep = '\t', quote = F, row.names = F, col.names = F)
write.table(de_genes_non.dune_noLFCthreshold$Ha412_gene,
            file = "analysis/GO_analysis/study_DE_genes_non-dune_noLFCthreshold.txt",
            sep = '\t', quote = F, row.names = F, col.names = F)

#### Diff expression manhattan plot ####
# Filter gff for only the mRNA, ncRNA, and rRNA (has same positions as 'gene';
# there are 8 ncRNAs in or DE gene list, notably)
transcripts <- filter(gff, type %in% c("mRNA", "ncRNA", "rRNA"))

# get the cumulative position of each gene and then
# the max gene end position on each chrom
chrom_lengths <- read.table("data/ref_genome_Ha412HO/chrom_sizes_Ha412HOv2.0.txt",
                            col.names = c("chrom", "length")) %>%
  mutate(length=as.numeric(length))

cum_lengths <- chrom_lengths %>%
  mutate(cumstart=cumsum(length)-length,
         cumend=cumsum(length))

# get inversion regions
inversion_regions <- read.table("analysis/inversions/inversion_regions.txt") %>%
  dplyr::select(V2,V3,V4,V9) %>%
  rename(chrom=V2, start=V3, end=V4, name=V9) %>%
  filter(name %in% c("pet05.01", "pet09.01", "pet11.01", "pet17.01")) %>%
  left_join(cum_lengths[c(1,3)]) %>%
  mutate(inv_cumstart=start+cumstart, inv_cumend=end+cumstart)

# expression manhattan plot with no LFC threshold on DE test
expr_manhattan_df_noLFCthreshold <- left_join(data.frame(de_results_Shrink) %>%
                                                tibble::rownames_to_column("Ha412_gene"),
                                              dplyr::select(transcripts,
                                                            Ha412_gene,
                                                            chrom,
                                                            start,
                                                            end,
                                                            strand)) %>%
  mutate(Ha412_gene=gsub("mRNA","gene",Ha412_gene)) %>%
  mutate(Ha412_gene=gsub("ncRNA","gene",Ha412_gene)) %>%
  inner_join(cum_lengths) %>% #inner join gets rid of non-chromosomal contigs
  mutate(gene_cumstart = start + cumstart, sig=if_else(padj<.05, "yes", "no"))

#for( i in 1:nrow(expr_manhattan_df_noLFCthreshold)){
#  if( !is.na(expr_manhattan_df_noLFCthreshold$padj[i]) &
#      expr_manhattan_df_noLFCthreshold$padj[i] == 0 ){
#    expr_manhattan_df_noLFCthreshold$padj[i] <- 6.996235e-204 #next smallest FDR
#  }
#}

sig_data_expr <- expr_manhattan_df_noLFCthreshold %>% 
  subset(padj < .05)
notsig_data_expr <- expr_manhattan_df_noLFCthreshold %>% 
  subset(padj >= .05) %>%
  group_by(chrom) %>% 
  sample_frac(0.75)

notsig_remove_expr <- subset(notsig_data_expr,abs(log2FoldChange) < .05) %>%
  sample_frac(0.95) #remove 95% of LFC <.05 non-significant genes.
notsig_keep_expr <- !(notsig_data_expr$Ha412_gene %in% notsig_remove_expr$Ha412_gene)
notsig_data_expr <- notsig_data_expr[notsig_keep_expr,]

expr_manhattan_df_noLFCthreshold_reduced <- bind_rows(sig_data_expr,
                                                      notsig_data_expr)

# get middle position of each chromosome for x axis label
axisdf_expr <- expr_manhattan_df_noLFCthreshold_reduced %>%
  group_by(chrom) %>%
  summarize(center=(max(gene_cumstart) + min(gene_cumstart)) / 2) %>%
  mutate(chrom=c("1","2","3","4","5",
                 "6","7","8","9","10",
                 "11","12","13","14","15",
                 "16","17"))

expr_manhattan_plot <- ggplot(expr_manhattan_df_noLFCthreshold_reduced,
                              aes(x=gene_cumstart, y=log2FoldChange)) +
                                  #alpha=as.factor(sig)=="yes")) +
                                  #color=-log10(padj))) +
  # alternate shading of chromosomes
  geom_rect(data=cum_lengths,
            aes(xmin=cumstart, xmax=cumend,
                ymin=min(expr_manhattan_df_noLFCthreshold_reduced$log2FoldChange),
                ymax=max(expr_manhattan_df_noLFCthreshold_reduced$log2FoldChange),
                fill=as.factor(chrom)),
            inherit.aes = F,
            alpha=0.5) +
  scale_fill_manual(values = rep(c("light grey", "white"), 17), guide = "none") +
  
  # mark inversion regions
  geom_rect(data=inversion_regions,
            aes(xmin=inv_cumstart, xmax=inv_cumend,
                ymin=min(expr_manhattan_df_noLFCthreshold_reduced$log2FoldChange),
                ymax=max(expr_manhattan_df_noLFCthreshold_reduced$log2FoldChange)),
            fill="red", alpha=0.25, inherit.aes = F) +

  # add each gene
  geom_point(aes(size=as.factor(sig=="yes")), shape=1) + #alpha=0.75, aes(size=-log10(padj))) +
  #scale_alpha_manual(name="FDR < .05",values = c(.1,1), labels=c("no", "yes")) +
  scale_size_manual(name="FDR < .05", values = c(1,4), labels=c("no", "yes")) +
  #scale_color_gradientn(colors=c("grey","turquoise","blue","navyblue","black")) +
  
  # custom axes
  scale_x_continuous(label = axisdf_expr$chrom,
                     breaks = axisdf_expr$center, expand=c(.01,.01)) +
  scale_y_continuous(expand=c(.01,.01)) +
  coord_cartesian(clip = 'off') +
  
  # custom theme
  theme_bw(base_size=18) +
  theme(legend.position = "none",
        #axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),
        #text = element_text(size=24),
        axis.line.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x=element_blank(),
        #plot.title = element_text(size=18),
        #axis.text.y = element_text(size=18),
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +

  labs(x="", y="log2FC")

expr_manhattan_plot
ggsave("figures/expr_manhattan.png",plot = expr_manhattan_plot, device = "png",
       width = 8, height = 6, units = "in", dpi = 300)

#### how many DE genes on each chromosome? ####
chrom_lengths <- read.table("data/ref_genome_Ha412HO/chrom_sizes_Ha412HOv2.0.txt",
                            col.names = c("chrom", "length"))
de_genes_per_chrom <- de_results %>% left_join(dplyr::select(transcripts,
                                                             ID, chrom)) %>%
  subset(padj<.05) %>%
  group_by(chrom) %>%
  filter(!grepl("Chr00", chrom)) %>%
  summarize(num_de_genes=n()) %>%
  right_join(chrom_lengths) %>%
  mutate(num_de_genes_norm = (num_de_genes/length)*1000000)
de_genes_per_chrom$chrom <- c("01","02","03","04","05","06","07","08","09",
                              "10","11","12","13","14","15","16","17")

DE_genes_per_chrom_plot <- ggplot(data=de_genes_per_chrom,
                                  
                                  aes(x=chrom, y=num_de_genes_norm, fill=chrom)) +
  geom_bar(stat = "identity", alpha=0.9) +
  theme_classic() +
  scale_fill_manual(values = rep(c("cornflowerblue", "grey50"), unique(length(DE_genes_per_chrom_plot)))) +
  labs(x="Chromosome", y="DE genes per 1Mb") +
  theme(text = element_text(size=36),
        legend.position = "none",
        legend.title = element_blank())
  
#### individual genes expression plots ####
norm_expr_df <- normalized_counts[,order(colnames(normalized_counts))] %>% 
  rownames_to_column("Ha412_gene") %>%
  mutate(Ha412_gene=gsub(".*RNA","gene",Ha412_gene)) %>%
  pivot_longer(2:25, names_to = "sample", values_to = "normalized_expression") %>%
  mutate(ecotype=if_else(grepl("non", sample), "Non-dune", "Dune"),
         Ha412_gene=gsub("gene:","",Ha412_gene))
  

# ATPase alpha subunit. ATCG00120. All 3 sunflower homologs are DE:
# "Ha412HOChr13g0618161","Ha412HOChr17g0827191", "Ha412HOChr17g0827211"
# Just plot the two more expressed homologs
ATPalpha_homologs <- c("Ha412HOChr13g0618161",
                       "Ha412HOChr17g0827191")

# ATCG00480 is chloroplast-encoded gene for beta subunit of ATP synthase.
# Sunflower appears to have multiple copies of this gene in the nuclear genome.
# All four DE but just two are mod to highly expressed.
# "Ha412HOChr01g0005241", "Ha412HOChr13g0616301",
# "Ha412HOChr08g0354771", "Ha412HOChr07g0317901"
ATPbeta_homologs <- c("Ha412HOChr01g0005241",
                "Ha412HOChr13g0616301")

# AT4G04640 i.e. ATPC1 encodes chloroplast ATP synthase subunit gamma
# Three homologs, none are significantly DE, but two show a trend of non-dune upreg
# "Ha412HOChr02g0069731" very lowly expressed.
ATPgamma_homologs <- c("Ha412HOChr05g0200861", 
                       "Ha412HOChr15g0702251")

# AT4G09650 encodes the chloroplast ATP synthase delta subunit. 
# two of the four homologs are DE: Ha412HOChr02g0070001 and Ha412HOChr05g0243771.
# "Ha412HOChr05g0243771", "Ha412HOChr10g0431791"
# "Ha412HOChr01g0001311", "Ha412HOChr02g0070001" are very lowly expressed though.
ATPdelta_homologs <- c("Ha412HOChr05g0243771",
                       "Ha412HOChr10g0431791")

# ATCG00470 subunit epsilon. Both are significantly DE but
# Ha412HOChr07g0308261 is very lowly expressed
ATPepsilon_homologs <- c("Ha412HOChr07g0308261", 
                         "Ha412HOChr13g0618521")

# ATCG00150 aka ATPI Encodes a subunit of ATPase complex CF0,
# which is a proton channel that supplies the proton motive force
# to drive ATP synthesis by CF1 portion of the complex.
# all 3 are significantly DE though Ha412HOChr13g0620231 is lowly expressed
ATPi_homologs <- c("Ha412HOChr12g0569601",
                   "Ha412HOChr13g0620221")

ATPase_list <- list(ATPalpha_homologs, ATPbeta_homologs, ATPgamma_homologs,
                    ATPdelta_homologs, ATPepsilon_homologs, ATPi_homologs)

ATPase_plot_list <- list()
i <- 1
for( subunit in ATPase_list ){
  plot <- ggplot(data=norm_expr_df %>% filter(Ha412_gene %in% subunit),
         aes(x=ecotype, y=normalized_expression, shape=ecotype, fill=ecotype)) +
    facet_wrap("Ha412_gene") +
    #geom_point(size=5) +
    #geom_violin(alpha=.25) +
    geom_point(size=3, alpha=.75, position = position_jitter(width = .2)) +
    geom_boxplot(alpha=0.25, outlier.shape=NA) +
    theme_bw(base_size = 12) +
    theme(legend.position = "none",
          strip.text = element_text(size = 7)) +
    scale_fill_manual(values = c("gold2", "forestgreen")) +
    scale_color_manual(values = c("gold2", "forestgreen")) +
    scale_shape_manual(values=c(24,21)) +
    labs(x="Ecotype", y="Normalized expression")
  ATPase_plot_list[[i]] <- plot
  i <- i + 1
}
ATPase_plot_list[[1]] <- ATPase_plot_list[[1]] +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) + 
  labs(x="",y="Norm. expression")#1x3

ATPase_plot_list[[2]] <- ATPase_plot_list[[2]] +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(x="", y="")

ATPase_plot_list[[3]] <- ATPase_plot_list[[3]] +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(x="",y="Norm. expression")

ATPase_plot_list[[4]] <- ATPase_plot_list[[4]] + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(x="", y="")

ATPase_plot_list[[5]] <- ATPase_plot_list[[5]] +
  labs(y="Norm. expression")

ATPase_plot_list[[6]] <- ATPase_plot_list[[6]] +
  labs(y="")

plot_ATPase_expr <- (ATPase_plot_list[[1]] | ATPase_plot_list[[2]]) /
  (ATPase_plot_list[[3]] | ATPase_plot_list[[4]]) /
  (ATPase_plot_list[[5]] | ATPase_plot_list[[6]])

ggsave("figures/plot_ATPase_expr_raw.v2.pdf", plot=plot_ATPase_expr, device = "pdf",
       height = 145.833, width = 175, dpi = 300, units = "mm")

# GLH17
plot_GLH17_expr <- ggplot(data=norm_expr_df %>%
                            filter(Ha412_gene %in% c("Ha412HOChr02g0088581")),
               aes(x=ecotype, y=normalized_expression,
                   shape=ecotype, fill=ecotype)) +
  geom_point(size=5, alpha=.75, position = position_jitter(width = .1)) +
  geom_boxplot(alpha=0.25, outlier.shape = NA) +
  theme_bw(base_size = 18) +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("gold2", "forestgreen")) +
  scale_color_manual(values = c("gold2", "forestgreen")) +
  scale_shape_manual(values=c(24,21)) +
  labs(x="Ecotype", y="Norm. expression")

# CESA6
plot_CESA6_expr <- ggplot(data=norm_expr_df %>%
                            filter(Ha412_gene %in% c("Ha412HOChr05g0243281")),
                          aes(x=ecotype, y=normalized_expression,
                              shape=ecotype, fill=ecotype)) +
  #facet_wrap("Ha412_gene") +
  geom_point(size=5, alpha=.75, position = position_jitter(width = .1)) +
  geom_boxplot(alpha=0.25, outlier.shape = NA) +
  theme_bw(base_size = 18) +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_fill_manual(values = c("gold2", "forestgreen")) +
  scale_color_manual(values = c("gold2", "forestgreen")) +
  scale_shape_manual(values=c(24,21)) +
  labs(x="Ecotype", y="Normalized expression")

# EMB2004; AT1G10510
plot_EMB2004_expr <- ggplot(data=norm_expr_df %>%
                            filter(Ha412_gene %in% c("Ha412HOChr17g0820321")),
                          aes(x=ecotype, y=normalized_expression,
                              shape=ecotype, fill=ecotype)) +
  #facet_wrap("Ha412_gene") +
  geom_point(size=5, alpha=.75, position = position_jitter(width = .1)) +
  geom_boxplot(alpha=0.25, outlier.shape = NA) +
  theme_bw(base_size = 18) +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_fill_manual(values = c("gold2", "forestgreen")) +
  scale_color_manual(values = c("gold2", "forestgreen")) +
  scale_shape_manual(values=c(24,21)) +
  labs(x="Ecotype", y="Normalized expression")

# SCPL13
plot_SCPL13_expr <- ggplot(data=norm_expr_df %>%
                             filter(Ha412_gene %in% c("Ha412HOChr05g0243281")),
                           aes(x=ecotype, y=normalized_expression,
                               shape=ecotype, fill=ecotype)) +
  facet_wrap("Ha412_gene") +
  geom_point(size=5, alpha=.75, position = position_jitter(width = .1)) +
  geom_boxplot(alpha=0.25, outlier.shape = NA) +
  theme_bw(base_size = 18) +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("gold2", "forestgreen")) +
  scale_color_manual(values = c("gold2", "forestgreen")) +
  scale_shape_manual(values=c(24,21)) +
  labs(x="Ecotype", y="Normalized expression") 

# SCPL13
plot_SCPL13_expr <- ggplot(data=norm_expr_df %>%
                             filter(Ha412_gene %in% c("Ha412HOChr15g0742211")),
                           aes(x=ecotype, y=normalized_expression,
                               shape=ecotype, fill=ecotype)) +
  facet_wrap("Ha412_gene") +
  #geom_point(size=5) +
  #geom_violin(alpha=.25) +
  geom_point(size=5, alpha=.75, position = position_jitter(width = .1)) +
  geom_boxplot(alpha=0.25, outlier.shape = NA) +
  theme_bw(base_size = 18) +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("gold2", "forestgreen")) +
  scale_color_manual(values = c("gold2", "forestgreen")) +
  scale_shape_manual(values=c(24,21)) +
  labs(x="Ecotype", y="Normalized expression") 

# ATO (SF3a60, AT5G06160)
# two Ha12HO homologs: Ha412HOChr09g0395181 and Ha412HOChr09g0395201; 
# the latter is strongly DE; the former is far more highly expressed, but not DE
plot_ATO_expr <- ggplot(data=norm_expr_df %>%
                             filter(Ha412_gene %in% c("Ha412HOChr09g0395201")),
                           aes(x=ecotype, y=normalized_expression,
                               shape=ecotype, fill=ecotype)) +
  facet_wrap("Ha412_gene") +
  geom_point(size=3, alpha=.75, position = position_jitter(width = .1)) +
  geom_boxplot(alpha=0.25, outlier.shape = NA) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none",
        strip.text = element_text(size = 6)) +
  scale_fill_manual(values = c("gold2", "forestgreen")) +
  scale_color_manual(values = c("gold2", "forestgreen")) +
  scale_shape_manual(values=c(24,21)) +
  labs(x="", y="Normalized expression")

# SUS2 (AT1G80070).
# Seven Ha412HO homologs, but just 3 are DE:
# Ha412HOChr02g0050341, Ha412HOChr16g0759321, Ha412HOChr16g0759311
# Note that the non-DE SUS2 homolog Ha412HOChr04g0149001 
# is the most highly expressed homolog at 7322 base mean;
# the DE homolog Ha412HOChr02g0050341 has next highest base mean of 213.
plot_SUS2_expr <- ggplot(data=norm_expr_df %>%
                          filter(Ha412_gene %in% c("Ha412HOChr02g0050341")),
                        aes(x=ecotype, y=normalized_expression,
                            shape=ecotype, fill=ecotype)) +
  facet_wrap("Ha412_gene") +
  geom_point(size=3, alpha=.75, position = position_jitter(width = .1)) +
  geom_boxplot(alpha=0.25, outlier.shape = NA) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none",
        strip.text = element_text(size = 6)) +
  scale_fill_manual(values = c("gold2", "forestgreen")) +
  scale_color_manual(values = c("gold2", "forestgreen")) +
  scale_shape_manual(values=c(24,21)) +
  labs(x="Ecotype", y="")

# CWC22 homolog splice factor
plot_CWC22_expr <- ggplot(data=norm_expr_df %>%
                            filter(Ha412_gene %in% c("Ha412HOChr08g0341371",
                                                     "Ha412HOChr08g0341411")),
                          aes(x=ecotype, y=normalized_expression,
                              shape=ecotype, fill=ecotype)) +
  facet_wrap("Ha412_gene") +
  geom_point(size=3, alpha=.75, position = position_jitter(width = .1)) +
  geom_boxplot(alpha=0.25, outlier.shape = NA) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none",
        strip.text = element_text(size = 6)) +
  scale_fill_manual(values = c("gold2", "forestgreen")) +
  scale_color_manual(values = c("gold2", "forestgreen")) +
  scale_shape_manual(values=c(24,21)) +
  labs(x="", y="")

plot_splice_factor_expr <- plot_ATO_expr + plot_SUS2_expr + plot_CWC22_expr +
  plot_layout(widths = c(1,1,1.75))
ggsave("figures/plot_splice_factor_expr_raw.pdf",
       plot = plot_splice_factor_expr,
       device ="pdf", height = 65.625, width = 175, units ="mm", dpi = 300)

#### Volcano plot ####
pval_outliers <- c("gene:Ha412HOChr01g0005241", "gene:Ha412HOChr13g0616301",
                     "gene:Ha412HOChr17g0815401", "gene:Ha412HOChr08g0354771")

data.frame(de_results_Shrink) %>%
  rownames_to_column("Ha412_gene") %>% 
  mutate(Ha412_gene=gsub(".*RNA","gene",Ha412_gene)) %>%
  filter(!Ha412_gene %in% pval_outliers) %>%
  # plot
  ggplot(aes(x=log2FoldChange, y=-log10(pvalue))) +
  geom_point(size=3, alpha=0.5) +
  geom_hline(yintercept = 2.1, color="red", lty=2) +
  theme_bw() +
  theme(text = element_text(size=24),
        axis.text = element_text(size=18))

#### MA plot ####
#MA_plot <- ggplot(data.frame(de_results_lfc1_Shrink), aes(x=baseMean, y=log2FoldChange, color=svalue<.005)) +
#  geom_point(alpha=0.5, size=3, show.legend = F) +
#  scale_color_manual(values = c('TRUE'='cornflowerblue', 'FALSE'='grey50')) +
#  scale_x_log10() +
#  theme_classic()
#
#DESeq2::plotMA(de_results_Shrink)
#DESeq2::plotMA(DE)
