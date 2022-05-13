# Differential Expression analysis using DESeq2
# http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
BiocManager::install("DESeq2")
BiocManager::install("apeglm")
library(DESeq2)
library(tidyr)
library(dplyr)
library(tibble)
library(ggplot2)
library(patchwork)
library(apeglm)

#?DESeqDataSetFromHTSeqCount

sample_table <- read.table("analysis/DESeq2/sample_table.txt", header=T)
sample_table$habitat <- as.factor(sample_table$habitat)
counts <- read.table("analysis/DESeq2/htseq-count_out/htseq-count_results.txt", row.names = 1)
names(counts) <- paste(sample_table$habitat,"_",sample_table$sample_no, sep = "")

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = sample_table,
                              design = ~ habitat)

# Filter out lowly expressed genes, where total read count is less than 10 across all samples combined
keep <- rowSums(counts(dds)) >= 10
low_expressed_genes <- data.frame(keep) %>% #make a df for filtering of the DEXSeq results in analyze_DEU_DEXSeq.R
  subset(keep=="FALSE") %>%
  rownames_to_column("groupID") %>%
  mutate(groupID=gsub("mRNA:","",groupID))
dds <- dds[keep,]
# re-level the habitat factor to make non_dune samples the reference (i.e. the denominator in log Fold change calculation, log2(dune/non-dune))
dds$habitat <- relevel(dds$habitat, ref = "non-dune")

# Analyze differential expression
DE <- DESeq(dds)
save(DE, file="analysis/DESeq2/DESeq2_results.Robj") #save so we don't have to repeat this step
load("analysis/DESeq2/DESeq2_results.Robj")

# Extract normalized counts. I don't think I end up using this for anything.
normalized_counts <- as.data.frame(DESeq2::counts(DE, normalized=TRUE))

# rlog transformation. for expression PCA. needs to be applied WITHOUT normalization bc it requires expression values to be integers
rlog_counts <- rlog(DESeq2::counts(DE, normalized=FALSE)) #34754 genes
# save to file for use in PCA
write.csv(rlog_counts, file = "analysis/pca/rlog-transformed_expression_counts.csv", quote = F, row.names = F)

#### DE results ####
# Provide contrast argument to set dune habitat in the numerator. Set LFC threshold of 1, which means a particular having twice as much expression in one ecotype vs the other. ALso, set significance threshold to .05 
de_results_p05 <- results(DE, contrast = c("habitat","dune","non-dune"), alpha = .05, cooksCutoff = F) #no LFC threshold
de_results_lfc1_p05 <- results(DE, contrast = c("habitat","dune","non-dune"), alpha = .05, lfcThreshold = 1, cooksCutoff = F) #threshold of LFC=1

de_results_p05 <- data.frame(de_results_p05) %>%
  rownames_to_column("id") %>%
  arrange(desc(abs(log2FoldChange))) %>%
  subset(padj<.05)

de_results_lfc1_p05 <- data.frame(de_results_lfc1_p05) %>%
  rownames_to_column("id") %>%
  arrange(desc(abs(log2FoldChange))) %>%
  subset(padj<.05)

# write list of genes to file
write.table(dplyr::select(de_results_p05, id), file = "analysis/DESeq2/DE_genes_padj.05.tsv", sep = '\t', quote = F, row.names = F)
write.table(dplyr::select(de_results_lfc1_p05, id), file = "analysis/DESeq2/DE_genes_LFC1_padj.05.tsv", sep = '\t', quote = F, row.names = F)


# shrink Log Fold Change for visualization and ranking? I think this makes visualization even more difficult bc a lot of our "good" DE genes above the abs(LFC)=1 threshold appear to have <1 LFC after the shrinkage. 
de_results_lfcS <- data.frame(lfcShrink(DE, coef=c("habitat_dune_vs_non_dune"), type="apeglm", lfcThreshold = 1)) %>%
  rownames_to_column("id") %>%
  arrange(desc(abs(log2FoldChange)))
#de_results_lfcS <- de_results_lfcS %>% left_join(dplyr::select(de_results, id, padj))

#### Map differential expression across the genome ####
# Read-in our annotations, for the gene positions. 
gff <- read.table("data/ref_genome_Ha412HO/HAN412_Eugene_curated_v1_1.gff3", col.names = c("chrom", "source", "type", "start", "end", "score", "strand", "phase", "attributes"))

# Filter for only the mRNA (has same positions as 'gene')
transcripts <- filter(gff, type=="mRNA") %>%
  # Parse the 'attribute' column. For now we will only keep
  # the first four fields of the attributes
  separate(attributes, into = c("name", "id", "parent", "GO_term"),
                           sep = ";") %>%
  # edit the 'ID' attribute to match the DE results file
  mutate(id=sub("ID=", "", id)) 

# get the cumulative position of each gene and then the max gene end position on each chrom
cumulative_expr <- transcripts %>% 
  filter(!grepl("Chr00", chrom)) %>%
  group_by(chrom) %>% 
  summarise(max_bp=max(c_across(start:end))) %>% 
  mutate(bp_add = lag(cumsum(as.numeric(max_bp)), default = 0)) %>% 
  dplyr::select(chrom, bp_add)

expr_manhattan_df <- left_join(de_results, dplyr::select(transcripts, id, chrom, start, end, strand)) %>%
  inner_join(cumulative_expr) %>% #inner join gets rid of non-chromosomal contigs. There were 10 sig DE genes on these unmapped contigs. 
  rowwise() %>%
  # get the cumulative start position of each gene
  # and yes/no if it is significantly differentially expressed
  mutate(bp_cumulative = start + bp_add, sig=if_else(padj<.05, "yes", "no"))

# split apart significant DE genes and nonsignifcant DE genes so we can downsample the nonsig data, for prettier plot
sig_data <- expr_manhattan_df %>% 
  subset(padj < 0.05)
notsig_data <- expr_manhattan_df %>% 
  subset(padj >= 0.05) %>%
  group_by(chrom) %>% 
  sample_frac(0.1)

expr_manhattan_df_reduced <- bind_rows(sig_data, notsig_data) %>%
  filter(!id %in% c("mRNA:Ha412HOChr16g0770771"))

# get middle position of each chromosome for axis label
axis_set_expr <- expr_manhattan_df %>%
  group_by(chrom) %>%
  summarize(max=max(bp_cumulative), min=min(bp_cumulative)) %>%
  mutate(center=(min+max)/2)
axis_set_expr$chrom <- c("01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17")

expr_manhattan_plot <- ggplot(expr_manhattan_df_reduced, aes(x=bp_cumulative, y=log2FoldChange, color=as.factor(chrom), alpha=as.factor(sig)=="yes")) +
  geom_point(size=3) +
  scale_x_continuous(label = axis_set_expr$chrom, breaks = axis_set_expr$center, expand=c(0.0175,0.0175)) +
  scale_color_manual(values = rep(c("cornflowerblue", "grey50"), unique(length(axis_set_expr$chrom)))) +
  scale_alpha_manual(values = c(.1,1)) +
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(size = 36)) +
  geom_hline(yintercept = -1, linetype = "dashed") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(x="")


#### how many DE genes on each chromosome? ####
chrom_lengths <- read.table("data/ref_genome_Ha412HO/chrom_sizes_Ha412HOv2.0.txt", col.names = c("chrom", "length"))
de_genes_per_chrom <- de_results %>% left_join(dplyr::select(transcripts, id, chrom)) %>%
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

FIG_2 <- plot_spacer() / Fst_manhattan_plot | (expr_manhattan_plot / DE_genes_per_chrom_plot)
ggsave("figures/FIG_2.jpg", plot=FIG_2, width=56, height = 28, dpi = 600, units = "cm", bg = "white")
  
#### MA plot ####
MA_plot <- ggplot(de_results, aes(x=baseMean, y=log2FoldChange)) +
  geom_point(aes(color=padj<0.05), pch=1, alpha=0.75, show.legend = F) +
  scale_color_manual(values = c('TRUE'='#FFC300', 'FALSE'='grey')) +
  scale_x_log10()

plotMA()
DESeq2::plotMA(DE)



