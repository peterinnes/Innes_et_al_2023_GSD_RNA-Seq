# script to analyze isoform expression/proportions and parents_diffv2 results 
library(dplyr)
library(ggplot2)
library(compositions)

#### read-in list of "good" alternatively spliced genes and isoforms,
# from the parents_diff_v2.py pipeline 
spliced_genes <- read.table("analysis/diff_iso_out/good_isoforms.post_consolidate_filtered.2022-10-10.txt",
                            col.names = c("gene","isoforms")) %>%
  separate(isoforms, into = paste0(1:10),
           sep = ",") %>%
  pivot_longer(cols = 2:11, names_to = "isoform_num", values_to = "isoform_id") %>%
  na.omit() %>%
  dplyr::select(gene,isoform_id)

head(spliced_genes)

# need to deal with TRINITY 'genes' and their isoforms that have been split into two 'genes'
# i.e. 'TRINITY...g1' becomes two genes 'TRINIRT...g1.0' and 'TRINITY...g1.1'

#### read-in and reformat TPM data

isoform_tpm_df <- read.table("data2/RSEM_out/isoform_tpm_matrix.txt",
                              header = T) %>%
  rename(transcript_id="isoform_id")

isoform_clusters <- read.table("analysis/diff_iso_out/consolidated_isoform_clusters_with_rownames.txt",
                               col.names=c("cluster","alleles")) %>%
  separate(alleles, into = paste0(1:5),
           sep = ",") %>%
  pivot_longer(cols = 2:6, names_to = "allele", values_to = "isoform_id") %>%
  na.omit() %>%
  inner_join(spliced_genes) %>%
  mutate(cluster=paste0(gene,"_",cluster)) %>%
  relocate(gene, cluster)

# sum TPM values for alleles within a cluster, for each sample
isoform_tpm_mat_clusters <- left_join(isoform_clusters, isoform_tpm_mat) %>%
  group_by(cluster, gene) %>%
  summarise_if(is.numeric, sum) %>%
  relocate(gene,cluster) %>%
  rename(cluster="isoform_id") #rename to match the no_clusters df

# subset the TPM matrix for isoforms (transcripts) that do not have 
# alternate alleles i.e. do not form a 'cluster' with other isoforms.
isoform_tpm_mat_no_clusters <- anti_join(isoform_tpm_mat, isoform_clusters) %>%
  inner_join(spliced_genes) %>%
  relocate(gene, isoform_id)

# combine the clusters and no_clusters isoform tpm dataframes into our final tpm df
diff_iso_tpm_df <- data.frame(rbind(isoform_tpm_mat_clusters,
                                    isoform_tpm_mat_no_clusters)) %>%
  arrange(gene)

# pivot longer for ggplotting
diff_iso_tpm_df_long <- pivot_longer(diff_iso_tpm_df, cols = 3:26,
                                     names_to = "sample", values_to = "TPM") %>%
  mutate(ecotype=gsub("_.*","", sample))

ggplot(data = subset(diff_iso_tpm_df_long, gene=="TRINITY_DN0_c0_g1.0"),
       aes(x=ecotype, y=TPM, color=isoform_id)) +
  geom_point()

#### convert TPM values to proportion of total TPM? ####
tpm_totals <- diff_iso_tpm_df %>% group_by(gene) %>%
  summarise_if(is.numeric, sum) %>%
  pivot_longer(cols = 2:25, names_to = "sample", values_to = "gene_total_TPM")


# read-in list of genes with just two isoforms. subset to these genes moving forward
# so that we can calculate Percent Spliced In more easily, and compare to rMATS results
Trinity_genes_with_two_isos <- read.table("analysis/diff_iso_out/genes_with_two_isoforms.txt",
                                          col.names = "gene")

diff_iso_tpm_two_isos_df <- filter(diff_iso_tpm_df_long, gene %in% Trinity_genes_with_two_isos$gene)

# use simple "1" or "2" do distinguish the alternative isoforms, rather than
# the actual isoform IDs. This will allow us to pivot_longer
diff_iso_prop_df <- left_join(diff_iso_tpm_two_isos_df, tpm_totals) %>%
  mutate(proportion_TPM = TPM / gene_total_TPM,
         isoform_num = rep(c(rep("iso_1", 24), rep("iso_2",24)),
                           length(unique(gene))),
         sample_num = rep(rep(seq(1:12), 2), length(unique(gene))*2))

diff_iso_ilt_df <- dplyr::select(diff_iso_prop_df, gene, ecotype, sample,
                                 isoform_num, proportion_TPM) %>%
  pivot_wider(names_from = isoform_num, values_from = proportion_TPM) %>% as.data.frame()

dune_ilt_df <- subset(diff_iso_ilt_df, ecotype=="dune") %>%
  na.omit()
non.dune_ilt_df <- subset(diff_iso_ilt_df, ecotype=="non.dune") %>% 
  na.omit()

# split data frames into list of dataframes so we can do the ilr transform
# for each gene
dune_ilt_list <- split( dune_ilt_df , f = dune_ilt_df$gene )
non.dune_ilt_list <- split( non.dune_ilt_df, f = non.dune_ilt_df$gene)

dune_res <- list()
i <- 1
for( gene_df in dune_ilt_list ){
  info <- gene_df[1:3]
  info$ilt <- as.vector(ilt(gene_df$iso_1, gene_df$iso_2))
  dune_res[[i]] <- info
  i <- i + 1
}

non.dune_res <- list()
i <- 1
for( gene_df in  non.dune_ilt_list ){
  info <- gene_df[1:3]
  info$ilt <- as.vector(ilt(gene_df$iso_1, gene_df$iso_2))
  non.dune_res[[i]] <- info
  i <- i + 1
}

diff_iso_splice_df <- bind_rows(dune_res, non.dune_res) %>%
  dplyr::select(!ecotype) %>%
  pivot_wider(names_from=sample, values_from = ilt) %>%
  dplyr::select(!gene)

#### PCA #### 
# of the isometric-log-ratio tranformed isoform proportions
# for just the genes with two alternative isoforms

# remove genes with > 40% missingness
diff_iso_splice_pca_df <- diff_iso_splice_df %>%
  mutate(missing_perc = rowMeans(is.na(diff_iso_splice_df))) %>%
  filter(missing_perc<.4) %>%
  dplyr::select(!missing_perc)

# set missing values to the average PSI of each event (each row)
for( i in 1:nrow(diff_iso_splice_pca_df) ){
  for( j in 1:ncol(diff_iso_splice_pca_df) ){
    if( is.na(diff_iso_splice_pca_df[i,j]) ){
      diff_iso_splice_pca_df[i,j] <- rowMeans(diff_iso_splice_pca_df[i,], na.rm = T)
    }
  }
}

diff_iso_splice_pca <- rda(t(diff_iso_splice_pca_df), scale = F, center=T)

diff_iso_splice_pve_df <- data.frame(summary(
  eigenvals(diff_iso_splice_pca)))[2,1:12] %>%
  pivot_longer(values_to = "PVE",
               names_to="PC",
               cols = 1:12) %>%
  mutate(PVE=100*PVE)

diff_iso_splice_pca.site_sc <- data.frame(scores(diff_iso_splice_pca, choices = 1:2,
                                              scaling = 1, display = "wa")) %>%
  tibble::rownames_to_column("sample.id") %>%
  mutate(habitat=gsub("_.*","", sample.id))

# get convex hulls for splicing pca plot
hull_diff_iso_d <- filter(diff_iso_splice_pca.site_sc, habitat=="dune") %>%
  dplyr::slice(chull(PC1,PC2))
hull_diff_iso_nd <- filter(diff_iso_splice_pca.site_sc, habitat=="non.dune") %>%
  dplyr::slice(chull(PC1,PC2))

# plotting
plot_diff_iso_splice_pca <-  ggplot(data=diff_iso_splice_pca.site_sc,
                                 aes(x = PC1, y = PC2, shape=habitat,
                                     fill=habitat)) +
  geom_point(size=5, alpha=.75) +
  geom_polygon(data=hull_diff_iso_d, alpha=0.5) +
  geom_polygon(data=hull_diff_iso_nd, alpha=0.5) +
  
  labs(x=paste0("PC1 ","(",round(diff_iso_splice_pve_df[1,2], 2),"%)"),
       y=paste0("PC2 ", "(",round(diff_iso_splice_pve_df[2,2],2),"%)")) +#,
  #title="Alternative splicing") +
  theme_bw() +
  theme(text = element_text(size=24),
        legend.position = "none",
        legend.title = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.title = element_text(size=18),
        axis.text = element_text(size=18)) +
  geom_hline(yintercept = 0, lty = 2, col = "grey50") +
  geom_vline(xintercept = 0, lty = 2, col = "grey50") +
  scale_color_manual(values=c("gold2", "forestgreen")) +
  scale_fill_manual(values=c("gold2", "forestgreen")) +
  scale_shape_manual(values=c(24,21))

#### trying just PSI instead of ILT ####
tmp <- dplyr::select(diff_iso_ilt_df, gene, sample, iso_1) %>% pivot_wider(names_from = sample, values_from = iso_1) %>%
  dplyr::select(!gene) %>%
  na.omit()

diff_iso_splice_pca <- rda(t(tmp), scale = F, center=T)
            