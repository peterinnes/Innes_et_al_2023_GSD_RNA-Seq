#install.packages("venneuler")
#install.packages("VennDiagram")
install.packages("eulerr")
library(eulerr)
library(ggplot2)
library(ggvenn)
#library(venneuler)
#library(VennDiagram)

set1 <- read.table("analysis/GO_analysis/study_DE_genes_noLFCthreshold.txt")
set2 <- read.table("analysis/GO_analysis/study_DS_rMATS_genes.txt")
set3 <- read.table("analysis/diff_iso_out/HAN412HOV2_all_hits.signigicant_parents_diff_iso_results.2022-10-10.txt")
length(base::intersect(set1$V1,set2$V1))

length(intersect(set2$V1,set3$V1)) #intersect b/w rMATS and parents_diff
write.table(intersect(set2$V1,set3$V1),
            file = "analysis/GO_analysis/study_DS_parents_diff_rMATS_intersect_genes.txt",
            quote = F, row.names = F, col.names = F)

#intersect b/w rMATS embryo dev genes and parents_diff
DS_rMATS_parents_diff_embryo_dev_genes <- intersect(set3$V1,DS_embryo_dev_genes)


length(intersect(set1$V1, DS_rMATS_parents_diff_embryo_dev_genes))


#set3 <- read.table("analysis/GO_analysis/study_DEU_genes_noLFCthreshold_padj.05.txt")
#set4 <- read.table("analysis/GO_analysis/study_DS_parents_diff_genes.txt")

# with eulerr package
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
                              fills = list(fill = c("#317EC2", "grey50", "red"),
                                           alpha=0.5)
)
png("figures/venn_dia_DE_DS_raw.png", width = 6, height = 4.5, units = "in",
    res = 300)#, bg = "transparent")
p
dev.off()
png("figures/venn_dia_DE_DS_blank.png", width = 6, height = 4.5, units = "in",
    res = 300, bg = "transparent")
d
dev.off()

#### hypergeometric test for geneset overlap significance ####
# x = num of genes in common between two groups.
# n = num of genes in group 1.
# D = num of genes in group 2.
# N = total genes,
# dhyper(x,n,(N-n),D)
dhyper(232,5080,(32308-5080),1038) #overlap of DE genes and rMATS DS genes

# representation factor = x / expected # of genes.
# Expected num of genes = (n * D) / N 
rf <- 232 / ((1038 * 5080)/32308) # = 1.42

#dhyper(324,5080,(32311-5080),956) #overlap of DE genes and DEXSeq DS genes
#dhyper(90, 5080,(32311-5080), 416) #overlap of DE genes and parents_diff genes
#dhyper(29, 5080, (32311-5080), 118) #overlap of DE genes DEXseq+rMATS DS genes 

#overlap_de_ds <- inner_join(set1, set2) %>%
#  #inner_join(set3) %>%
#  #inner_join(set4) %>%
#  rename("Ha412_gene"=V1) %>%
#  left_join(Ha412_Ath_mappings)
#
#
#rownames(set1) <- set1$V1
#set1 <- rownames(set1)
#rownames(set2) <- set2$V1
#set2 <- rownames(set2)
##rownames(set3) <- set3$V1
##set3 <- rownames(set3)
##rownames(set4) <- set4$V1
##set4 <- rownames(set4)
##
#geneset <- list(DE=set1, DS=set2)
#
#
#venn_dia_DE_DS <- ggvenn(data=geneset, columns = c("DE", "DS"),
#                         show_percentage = F, set_name_size = 10, text_size = 12,
#                         fill_color = c("black", "grey50"))
#
#ggsave(filename = "figures/venn_dia_DE_DS.png", device = "png",
#       plot=venn_dia_DE_DS, height=6, width=9, dpi = 300,
#       units = "in", bg = "white")
