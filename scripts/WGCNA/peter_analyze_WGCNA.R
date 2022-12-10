library(DESeq2)
library(WGCNA)
library(dplyr)

#### read-in, filter, and transform raw count data ####
sample_table_deseq2 <- read.table("analysis/DESeq2/sample_table.txt", header=T)
sample_table_deseq2$habitat <- as.factor(sample_table_deseq2$habitat)
counts <- read.table("analysis/DESeq2/htseq-count_out/htseq-count_results.2022-6-26.txt",
                     row.names = 1)

names(counts) <- paste(sample_table_deseq2$habitat,"_",
                       sample_table_deseq2$sample_no, sep = "")

# filter the raw counts 
keep_wgcna <- which((rowMeans(as.matrix(counts)) >= 1)) #same prefilter as for DESeq2
#keep_wgcna <- which((rowMedians(as.matrix(counts)) >= 1))
#keep_wgcna <- which((rowSums(as.matrix(counts)) >= 24))

filtered_counts_wgcna <- counts[keep_wgcna,] %>%
  head(-3) #remove the last three rows, "__no_feature", "__ambiguous" "__alignment_not_unique"
  
dim(filtered_counts_wgcna)

rlog_filtered_counts_wgcna <- rlog(as.matrix(filtered_counts_wgcna)) %>%
  data.frame()

write.table(rlog_filtered_counts_wgcna,
            file="analysis/WGCNA/rlog_filtered_counts_wgcna.tsv",
            quote = F, row.names = T, sep = "\t")

#### reformat expression data into a list of lists for WGCNA ####
non_dune_expr_data <- data.frame(t(rlog_filtered_counts_wgcna %>%
  dplyr::select(starts_with("non.dune"))))
dune_expr_data <- data.frame(t(rlog_filtered_counts_wgcna %>% 
  dplyr::select(starts_with("dune"))))

multi_expr <- list(non_dune = list(data = non_dune_expr_data),
                   dune = list(data = dune_expr_data))
                   
checkSets(multi_expr) #check that formatting is correct

#### build networks ####
nd_network <- blockwiseModules(multi_expr[[1]]$data, maxBlockSize = 35000,
                              power = 16, minModuleSize = 30,
                              networkType = "signed", deepSplit = 2,
                              TOMType = "signed",mergeCutHeight = 0.25,
                              numericLabels = TRUE, minKMEtoStay = 0,
                              saveTOMs = FALSE,  verbose = 4, nThreads = 12)

d_network <- blockwiseModules(multi_expr[[2]]$data, maxBlockSize = 35000,
                              power = 16, minModuleSize = 30,
                              networkType = "signed", deepSplit = 2,
                              TOMType = "signed",mergeCutHeight = 0.25,
                              numericLabels = TRUE, minKMEtoStay = 0,
                              saveTOMs = FALSE, verbose = 4, nThreads = 12)

# module labels as colors for nondune set
nd_colors = labels2colors(nd_network$colors)
# show number and size of modules
table(nd_network$colors)

# module labels as colors for dune set
d_colors = labels2colors(d_network$colors)
# matching labels to those used in nondune set 
d_colors = matchLabels(d_colors, nd_colors)
# show number and size of modules
table(dune_network$colors)

#average module size
mean(table(nondune_network$colors))
mean(table(dune_network$colors))

#save results 
save(nondune_network, nondune_colors, dune_network, dune_colors, file = "data/Rdata/independent_networks.Rdata")


save(nd_network, file = "non-dune_wgcna_network.Rdata")

