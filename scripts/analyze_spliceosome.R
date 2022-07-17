library(lme4)
library(gtools)
spliceosome_model_df <- data.frame(sample.id=splice_pca.site_sc$sample.id,
                                   splice_PC1=splice_pca.site_sc$PC1) %>%
  inner_join(dplyr::select(snp_pca_df, sample.id, snp_PC1=EV1))
  

Ha412_spliceosome <- read.table("analysis/BLAST_out/uniq_top_hits_spliceosome_components_vs_HAN412.txt", col.names = "id")
Ha412_spliceosome <- rbind(Ha412_spliceosome, c("mRNA:Ha412HOChr09g0395201")) #adding in the splicing factor SF3a60 homolog that is a top DE gene

spliceosome_expression <- subset(normalized_counts, rownames(normalized_counts) %in% Ha412_spliceosome$id)

sample_names <- data.frame(sample_name=mixedsort(names(spliceosome_expression)), sample.id=mixedsort(spliceosome_model_df$sample.id, decreasing = T))

spliceosome_model_df <- sample_names %>% left_join(spliceosome_model_df)
rownames(spliceosome_model_df) <- spliceosome_model_df$sample_name
spliceosome_model_df <- spliceosome_model_df[-c(1,2)]

# merge expression data and PC data. merge by rownames
spliceosome_model_df <- merge(spliceosome_model_df, t(spliceosome_expression), by=0)

# get rid of the "mRNA:" prefix
names(spliceosome_model_df) <- gsub("mRNA:","", names(spliceosome_model_df))
head(spliceosome_model_df)

# pivot longer
#test <- spliceosome_model_df %>%
#  pivot_longer(cols = 4:127,names_to = "gene",values_to = "normalized_counts")

#### this model doesn't work but was trying to answer the question of whether expression of splicing genes explained the pattern of diferential splicing ####
cool_model <- lm(splice_PC1 ~ Ha412HOChr09g0395201 + snp_PC1, data = spliceosome_model_df)
summary(cool_model)

#### which spliceosomal genes are differentially expressed and/or differentially spliced? ####

# DE
Ha412_spliceosome$id %in% de_results_lfc1_s005_Shrink$id #none of them

# DS (DEXSeq)
Ha412_spliceosome$id %in% deu_genes$groupID # none of them?

# DS (rMATS)
Ha412_spliceosome$id %in% gsub("gene","mRNA",sig_AS_events$GeneID) #one of them! Ha412HOChr17g0855721.
subset(all_AS_events_deltaPSI, GeneID=="gene:Ha412HOChr17g0855721") #the significant DS event is A3SS
