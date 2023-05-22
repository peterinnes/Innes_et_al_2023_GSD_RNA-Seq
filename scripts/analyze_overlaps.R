library(ggplot2)
library(eulerr)

set1 <- read.table("analysis/GO_analysis/study_DE_genes_noLFCthreshold.txt")
set2 <- read.table("analysis/GO_analysis/study_DS_rMATS_genes.txt")
set3 <- read.table("analysis/diff_iso_out/2023-02-01/HAN412HOv2_all_hits_bit100.signigicant_parents_diff_iso_results.2023-02-05.txt")
set4 <- read.table("analysis/diff_iso_out/2023-02-01/HAN412HOv2_all_hits_bit100.AS_genes.2023-02-05.txt")

# Trinity DS genes matched to Ha412HOv2 (top hits) 
set5 <- read.table("analysis/diff_iso_out/2023-02-01/HAN412HOv2_top_hits.significant_parents_diff_iso_results.2023-02-05.txt") %>%
  dplyr::select(V2) %>%
  mutate(V2=gsub("mRNA","gene",V2))

# Trinity DS genes matched to HA412HOv2 (reciprocal top hits) 
set6 <- read.table("analysis/diff_iso_out/2023-02-01/HAN412HOv2_reciprocal_best_hits.significant_parents_diff_iso_results.2023-02-05.txt") %>%
  dplyr::select(V1) %>%
  mutate(V1=gsub("mRNA","gene",V1))

length(base::intersect(set1$V1,set2$V1))

# intersect b/w rMATS and diff_iso DS genes w/ relaxed BLAST (bit score > 100) is 289.
length(intersect(set2$V1,set3$V1))

# intersect b/w rMATS and diff_iso DS genes w/ "high confidence" BLAST (single top Ha412HO hit for each Trinity gene) is 100
length(intersect(set2$V1,set5$V2)) 

#intersect b/w SMith et al AS genes and rMATS AS genes is 4582
length(intersect(set4$V1, unique(all_AS_events$GeneID)))

write.table(intersect(set2$V1,set3$V1),
            file = "analysis/GO_analysis/study_DS_parents_diff_rMATS_intersect_genes.txt",
            quote = F, row.names = F, col.names = F)

# intersect b/w rMATS embryo dev genes and parents_diff
DS_embryo_dev_genes <- read.table("analysis/GO_analysis/results_DS_rMATS_GO_enrichment.txt",
                                  header=T, sep="\t") %>%
  filter(name=="embryo development ending in seed dormancy") %>%
  dplyr::select(study_items) 
DS_embryo_dev_genes <- unlist(strsplit(x = DS_embryo_dev_genes$study_items,split = ", "))

DS_rMATS_parents_diff_embryo_dev_genes <- intersect(set3$V1,DS_embryo_dev_genes)

length(intersect(set1$V1, DS_rMATS_parents_diff_embryo_dev_genes))

# with eulerr package. total DE = 5103; total DS = 1038; intersect = 232
venn_dia_DE_DS <- euler(c(DE=4871, DS=806, "DE&DS"=232))
p <- eulerr:::plot.euler(venn_dia_DE_DS,
                         labels = NA, #list(fontsize=36),
                         quantities = list(fontsize=30),
                         fills = list(fill = c("#317EC2", "red", "purple"),
                                      alpha=0.5)
)

d <- eulerr:::plot.euler(venn_dia_DE_DS,
                              labels = NA,
                              quantities = NA,
                              fills = list(fill = c("#317EC2", "red", "purple"),
                                           alpha=0.5)
)
png("figures/venn_dia_DE_DS_raw.png", width = 6, height = 4.5, units = "in",
    res = 300)#, bg = "transparent")
p
dev.off()
png("figures/venn_dia_DE_DS_blank.png", width = 6, height = 4.5, units = "in",
    res = 300)#, bg = "transparent")
d
dev.off()

# Euler diagrams for rMATS and Smith et al approach.
# rMATS DS: 1038 - (overlap = 289) = 749
# diff_iso (Smith et al) DS: 1299 - 289 = 1010
venn_dia_rMATS_diff_iso <- euler(c(rMATS=749, diff_iso=1010,
                                   "rMATS&diff_iso"=289))

e <- eulerr:::plot.euler(venn_dia_rMATS_diff_iso,
                    labels = NA,
                    quantities = NA,
                    fills = list(fill = c("red", "pink1", "hotpink3"),
                                 alpha=0.5))

png("figures/venn_dia_rMATS_diff_iso_blank.png", width = 6, height = 4.5,
    units = "in", res = 300)#, bg = "transparent")
e
dev.off()


# overlap between rMATS AS genes and Smith et al ('diff_iso') AS genes
venn_dia_rMATS_diff_iso_all_AS <- euler(c(rMATS=2044, diff_iso=1425,
                                          "rMATS&diff_iso"=4582))
f <- eulerr:::plot.euler(venn_dia_rMATS_diff_iso_all_AS,
                         quantities = NA,
                         labels = NA,
                         fills = list(fill = c("grey75", "black", "grey50"),
                                      alpha=.5))
pdf("figures/venn_dia_rMATS_diff_iso_all_AS_blank.pdf", width = 6, height = 4.5)
f
dev.off()



#### hypergeometric test for geneset overlap significance ####
# x = num of genes in common between two groups.
# n = num of genes in group 1.
# D = num of genes in group 2.
# N = total genes,
# dhyper(x,n,(N-n),D)
dhyper(232,5103,(32308-5103),1038) #overlap of DE genes and rMATS DS genes

# representation factor = x / expected # of genes.
# Expected num of genes = (n * D) / N 
rf <- 232 / ((1038 * 5103)/32308) # = 1.42

#### make table for the 21 DE–DS (high confidence) genes ####
# these are 21 genes that are DE and DS (rMATS) and also DS according to
# Smith et al approach, using only the reciprocal top blastn hit b/w 
# the transcriptome and the genome. 

de_ds_ds_genes <- intersect(set1$V1, intersect(set2$V1, set5$V2))

de_ds_ds_genes_results <- de_results_Shrink_df %>%
  filter(Ha412HOv2_gene %in% de_ds_ds_genes) %>%
  left_join(filter(all_AS_events_deltaPSI, FDR <.05) %>%
              rename(Ha412HOv2_gene=GeneID)) %>%
  dplyr::select(chrom, gene=Ha412HOv2_gene, Araport11_gene, baseMean, log2FoldChange,
                AS_event=ID, IncLevelDifference) %>%
  group_by(gene) %>%
  slice_max(abs(IncLevelDifference)) %>%
  arrange(desc(abs(IncLevelDifference))) %>%
  as.data.frame()

write.table(x = de_ds_ds_genes_results, file = "analysis/de_ds_ds_genes_results.tsv",
            row.names = F, quote = F, sep = '\t')

# DE and rMATS DS, plus DS according to Smith et al approach using all 
# blastn hits > 100
de_ds_ds_genes_low <- intersect(set1$V1, intersect(set2$V1, set3$V1))

#### look through the 53 "very high confidence" DS–DS genes ####
ds_ds_genes_highest <- data.frame(GeneID=intersect(set2$V1,set6$V1),
                                  DS_confidence = "highest")

ds_ds_genes_results <- filter(all_AS_events_deltaPSI, 
                              GeneID %in% ds_ds_genes) %>%
  filter(FDR<.05) %>%
  group_by(GeneID) %>%
  slice_max(abs(IncLevelDifference)) %>%
  #arrange(FDR, desc(abs(IncLevelDifference)))
  arrange(desc(abs(IncLevelDifference)))

head(ds_ds_genes_results, 20)

ds_ds_genes_high <- data.frame(GeneID=intersect(set2$V1,set3$V1),
                               DS_confidence=NA) %>%
  mutate(DS_confidence=if_else(GeneID %in% ds_ds_genes_highest$GeneID,
          "highest","high"))

ds_results_df <- all_AS_events_deltaPSI %>% 
  left_join(ds_ds_genes_high) %>%
  rename(Ha412HOv2_gene=GeneID, Araport11_gene=Ath_gene)


for (i in 1:length(ds_results_df$DS_confidence)) {
  if (is.na(ds_results_df$DS_confidence[i]) & 
      ds_results_df$Ha412HOv2_gene[i] %in% set2$V1) {
    ds_results_df$DS_confidence[i] <- 'rMATS only'
  }
  if (ds_results_df$FDR[i] >= .05) {
    ds_results_df$DS_confidence[i] <- NA
  }
}

# Extended Data S2
write.table(ds_results_df,
            file="analysis/DS_results.tsv", sep = '\t',
            quote = F, row.names = F)
