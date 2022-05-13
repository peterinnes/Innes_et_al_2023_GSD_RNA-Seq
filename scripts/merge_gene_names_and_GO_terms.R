#### parsing BLAST hits for GO analysis ####
library(dplyr)

# read-in file that has BLASTx matches between Ha412 genes and Araport11 genes
Ha412_Ath_mappings <- read.table("analysis/BLAST_out/top_hits_HAN412_vs_Araport11.txt", col.names = c("Ha412_gene", "Ath_gene"))

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
