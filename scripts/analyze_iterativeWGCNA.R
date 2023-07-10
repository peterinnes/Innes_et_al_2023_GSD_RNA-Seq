# analyze iterativeWGCNA results

#### get the expression (count) data ####
sample_table_deseq2 <- read.table("analysis/DESeq2/sample_table.txt", header=T)
sample_table_deseq2$habitat <- as.factor(sample_table_deseq2$habitat)
counts <- read.table("analysis/DESeq2/htseq-count_out/htseq-count_results.2022-6-26.txt",
                     row.names = 1)
names(counts) <- paste(sample_table_deseq2$habitat,"_",
                       sample_table_deseq2$sample_no, sep = "")
keep_vstringent_wgcna <- which(rowMeans(as.matrix(counts)) >= 10)
filtered_vstringent_counts_wgcna <- counts[keep_vstringent_wgcna,] %>%
  head(-3) %>% 
  mutate(d_n_zero=rowSums(dplyr::select(.,starts_with("dune")) == 0),
         nd_n_zero=rowSums(dplyr::select(.,starts_with("non-dune")) == 0)) %>%
  filter(d_n_zero <= 6 & nd_n_zero <= 6) %>%
  dplyr::select(!c(d_n_zero, nd_n_zero))rlog_filtered_vstringent_counts_wgcna <- rlog(as.matrix(filtered_vstringent_counts_wgcna)) %>%
  data.frame()
non_dune_expr_data <- data.frame(t(rlog_filtered_vstringent_counts_wgcna %>%
                                     dplyr::select(starts_with("non.dune"))))

dune_expr_data <- data.frame(t(rlog_filtered_vstringent_counts_wgcna %>% 
                                 dplyr::select(starts_with("dune"))))

num_sets <- 2 #two sets, dune and non dune
set_labels <- c("Non-dune", "Dune")
multi_iter_expr <- list(non_dune = list(data = non_dune_expr_data),
                   dune = list(data = dune_expr_data))

#### read-in results from iterativeWGCNA ####
non_dune_iter_res <- read.table("analysis/WGCNA/iterativeWGCNA/non_dune/final-membership.txt",
                                header = T)
dune_iter_res <- read.table("analysis/WGCNA/iterativeWGCNA/dune/final-membership.txt",
                            header = T)


#### assign module labels to colors ####
nd_iter_mods <- non_dune_iter_res$Module
names(nd_iter_mods) <- non_dune_iter_res$Gene
nd_iter_colors <- labels2colors(nd_iter_mods)

d_iter_mods <- dune_iter_res$Module
names(d_iter_mods) <- dune_iter_res$Gene
d_iter_colors <- labels2colors(d_iter_mods)

# match module labels b/w networks (ecotypes)
d_iter_colors <- matchLabels(d_iter_colors, nd_iter_colors)

# make a dataframe of the module colors and sizes
iter_mod_sizes_df <- data.frame(table(d_iter_colors)) %>%
  dplyr::rename(iter_mod_color=d_iter_colors) %>%
  mutate(Ecotype="Dune") %>%
  full_join(data.frame(table(nd_iter_colors)) %>% 
              dplyr::rename(iter_mod_color=nd_iter_colors) %>%
              mutate(Ecotype="Non-dune")) %>%
  dplyr::rename(iter_mod_size=Freq)

# linear model, mod size ~ Ecotype. no sig difference
fit_iter_mod_sizes <- lm(log(iter_mod_size) ~ Ecotype - 1, data=iter_mod_sizes_df)

# number of modules with the same labels 
length(intersect(unique(d_colors), unique(nd_colors)))

# function for unique modules
outersect <- function(x, y) {
  sort(c(setdiff(x, y),
         setdiff(y, x)))
}

outersect(unique(d_iter_colors), unique(nd_iter_colors))

# Determining unique modules to each network: 
# modules unique to non-dune
nd_iter_unique_mods = setdiff(unique(nd_iter_colors), unique(d_iter_colors))

# modules unique to dune 
d_iter_unique_mods = setdiff(unique(d_iter_colors), unique(nd_iter_colors))

# calculate module preservation
multi_iter_color = list("Non-dune" = nd_iter_colors, "Dune" = d_iter_colors)
iter_mod_preservation = modulePreservation(multi_iter_expr, multi_iter_color,
                                         dataIsExpr = TRUE, 
                                         networkType = "signed",
                                         referenceNetworks = 1,
                                         nPermutations = 100, randomSeed = 1,
                                         quickCor = 0, verbose = 3,
                                         maxModuleSize = 2000)

save(iter_mod_preservation, file = "data2/Rdata/iterWGCNA_module_preservation.vanilla_defaults.Rdata")

#### dendrogram ####
#TOMsimilarityFromExpr(dune_expr_data, networkType="signed", power = 18)