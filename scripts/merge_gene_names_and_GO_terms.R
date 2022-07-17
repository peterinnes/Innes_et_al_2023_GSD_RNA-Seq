#### parsing BLAST hits for GO analysis ####
library(dplyr)

# read-in file that has BLASTx matches between Ha412 genes and Araport11 genes
Ha412_Ath_mappings <- read.table("analysis/BLAST_out/top_hits_HAN412_vs_Araport11.txt", col.names = c("Ha412_gene", "Ath_gene"))
#Ha412_Ath_mappings <- read.table("analysis/BLAST_out/temp", col.names = c("Ha412_gene", "Ath_gene"))
Ha412_Ath_mappings$Ha412_gene <- gsub("mRNA","gene", Ha412_Ath_mappings$Ha412_gene)

#stringtie_Ath_mappings <- read.table("analysis/BLAST_out/top_hits_stringtie_vs_Araport11.txt", col.names = c("stringtie_gene", "Ath_gene"))

#stringtie_Ha412_mappings <- read.table("data/stringtie/merged_transcripts.stringtie.Ha412_gene_names.txt", col.names = c("stringtie_gene", "Ha412_gene"))

# read-in our reformatted GO associations file
Ath_GO_associations <- read.table("data/ref_genome_Araport11/ATH_GO_GOSLIM.id2gos.txt", sep = '\t', col.names = c("Ath_gene", "GO_terms"))
  
#test <- cSplit(Ath_GO_associations, "GO_terms", sep=";")

# merge by Ath gene name to get GO terms for Ha412 genes.
Ha412_GO_associations <- left_join(Ha412_Ath_mappings, Ath_GO_associations, by="Ath_gene") %>%
  dplyr::select(Ha412_gene, GO_terms)

## get rid of the prefix "mRNA" so it's just the gene id.
#Ha412_GO_associations$Ha412_gene <- gsub("mRNA:", "", Ha412_GO_associations$Ha412_gene)

head(Ha412_GO_associations)
dim(Ha412_GO_associations)

# write to file
write.table(Ha412_GO_associations, file = "analysis/GO_analysis/Ha412_GO_GOSLIM.id2gos.txt", quote = F, sep = "\t", row.names = F, col.names = F)

#### cluster GO terms with Revigo web app ####
# Revigo results table of DE genes upregulated in Dune ecotype
revigo_DE_dune <- read.table("analysis/GO_analysis/Revigo_BP_table.DE_dune.tsv", header = T)

# Revigo results table of all DS genes (all unique rMATS + DEXSeq significant DS/DEU genes)
revigo_DS_rMATS_DEU_union <- read.table("analysis/GO_analysis/Revigo_BP_table.DS_rMATS_DEU_union.tsv", header=T, na.strings = "null")


plot_revigo_DE_dune <-  ggplot(data=revigo_DE_dune, aes(x = PC_1, y = PC_0, size=Value, color=Value)) +
  scale_size_continuous(range = c(2,8)) +
  geom_point(alpha=.5) +
  geom_text_repel(label=ifelse(revigo_DE_dune$Value>10,as.character(revigo_DE_dune$Name),''), force = 5, force_pull = 0.5, direction = "both") +
  labs(x="MDS2 (semantic space x)", y="MDS1 (semantic space y)",
       title="DE (dune)") +
  theme_bw() +
  geom_hline(yintercept = 0, linetype="dashed", color="gray50", alpha=.5) +
  geom_vline(xintercept = 0, linetype="dashed", color="gray50", alpha=.5) +
  theme(text = element_text(size=36),
        legend.position = "none",
        legend.title = element_blank())
        #legend.position = "bottom",
        #legend.margin=margin(0,0,0,0),
        #legend.box.margin=margin(0,0,0,0))
plot_revigo_DE_dune
ggsave(filename = "figures/revigo_mds_DE_dune.png", device = "png", plot=plot_revigo_DE_dune, height=6, width=8, dpi = 600, units = "in", bg = "white")

plot_revigo_DS <-  ggplot(data=revigo_DS_rMATS_DEU_union, aes(x = PC_1, y = PC_0, size=Value, color=Value)) +
  scale_size_continuous(range = c(2,8)) +
  geom_point(alpha=.5) +
  geom_text_repel(label=ifelse(revigo_DS_rMATS_DEU_union$Value>5,as.character(revigo_DS_rMATS_DEU_union$Name),''), force = 5, force_pull = 0.5, direction = "both") +
  labs(x="MDS2 (semantic space x)", y="",
       title="DS") +
  theme_bw() +
  geom_hline(yintercept = 0, linetype="dashed", color="gray50", alpha=.5) +
  geom_vline(xintercept = 0, linetype="dashed", color="gray50", alpha=.5) +
  theme(text = element_text(size=36),
        legend.position = "none",
        legend.title = element_blank())
#legend.position = "bottom",
#legend.margin=margin(0,0,0,0),
#legend.box.margin=margin(0,0,0,0))
ggsave(filename = "figures/revigo_mds_DS.png", device = "png", plot=plot_revigo_DS, height=6, width=8, dpi = 600, units = "in", bg = "white")

 
