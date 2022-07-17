library(ggvenn)
library(ggplot2)

set1 <- read.table("analysis/GO_analysis/study_DE_genes_noLFCthreshold.txt")
set2 <- read.table("analysis/GO_analysis/study_DEU_genes_noLFCthreshold_padj.05.txt")
set3 <- read.table("analysis/GO_analysis/study_DS_rMATS_genes.txt")
set4 <- read.table("analysis/GO_analysis/study_DS_parents_diff_genes.txt")

overlap_de_ds <- inner_join(set1, set2) %>%
  inner_join(set3) %>%
  inner_join(set4) %>%
  rename(V1="Ha412_gene") %>%
  left_join(Ha412_Ath_mappings)


rownames(set1) <- set1$V1
set1 <- rownames(set1)
rownames(set2) <- set2$V1
set2 <- rownames(set2)
rownames(set3) <- set3$V1
set3 <- rownames(set3)
rownames(set4) <- set4$V1
set4 <- rownames(set4)

geneset <- list(DESeq2=set1, DEXSeq=set2, rMATS=set3, parents_diff=set4)

#venn_dia_DE_DS <- ggvenn(geneset, c("DEXSeq", "rMATS", "parents_diff", "DESeq2"), show_percentage = F, set_name_size = 10, text_size = 14, fill_color = c("#B9DBF4", "darkblue", "cornflowerblue", "grey50"))
venn_dia_DE_DS <- ggvenn(geneset, c("DEXSeq", "rMATS", "DESeq2"), show_percentage = F, set_name_size = 10, text_size = 14, fill_color = c("#B9DBF4", "cornflowerblue", "grey50"))

ggsave(filename = "figures/poster_venn_dia_DE_DS.png", device = "png", plot=venn_dia_DE_DS, height=8, width=8, dpi = 600, units = "in", bg = "white")

#### hypergeometric test for geneset overlap significance ####
dhyper(324,5080,(32311-5080),956) #overlap of DE genes and DEXSeq DS genes
dhyper(232,5080,(32311-5080),1038) #overlap of DE genes and rMATS DS genes
dhyper(90, 5080,(32311-5080), 416) #overlap of DE genes and parents_diff genes

dhyper(29, 5080, (32311-5080), 118) #overlap of DE genes DEXseq+rMATS DS genes 
