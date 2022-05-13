library(ggvenn)
library(ggplot2)

set1 <- as.vector(read.table("analysis/GO_analysis/study_genes_DE_LFC1_padj.05.tsv"))
rownames(set1) <- set1$V1
set1 <- rownames(set1)

set2 <- read.table("analysis/GO_analysis/study_genes_DEU_LFC1_padj.05.tsv")
rownames(set2) <- set2$V1
set2 <- rownames(set2)
geneset <- list(DE=set1, DS=set2)

venn_dia_DE_DEU <- ggplot() + 
  geom_venn(geneset, textsize = 10) +
  theme_void() +
  theme(legend.position = "none") +
  scale_fill_manual(values=c("grey50", "grey50"))

ggsave("figures/venn_dia_DE_DEU.jpg", plot=venn_dia_DE_DEU, width=17, dpi = 600, units = "cm", bg = "white")

