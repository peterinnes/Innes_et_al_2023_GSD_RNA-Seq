# script to analyze isoform expression/proportions and parents_diffv2 results 
library(dplyr)
library(ggplot2)
library(ggridges)
library(forcats)
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
  rename(isoform_id=transcript_id)

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
  pivot_longer(cols = 2:25, names_to = "sample",
               values_to = "gene_total_TPM") %>%
  mutate(gene_total_TPM = gene_total_TPM + 0.000001) #add small value to avoid zeros

# proportions df with all genes
diff_iso_all_prop_df <- left_join(diff_iso_tpm_df_long, tpm_totals) %>%
  mutate(proportion_TPM = TPM / gene_total_TPM,
         ecotype = recode(ecotype, "dune" = "Dune", "non.dune" = "Non-dune"))

# read-in list of genes with just two isoforms. subset to these genes moving forward
# so that we can calculate Percent Spliced In more easily, and more easily compare to rMATS results
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

diff_iso_ilr_df <- dplyr::select(diff_iso_prop_df, gene, ecotype, sample,
                                 isoform_num, proportion_TPM) %>%
  pivot_wider(names_from = isoform_num, values_from = proportion_TPM) %>% as.data.frame()

dune_ilr_df <- subset(diff_iso_ilr_df, ecotype=="dune") %>%
  na.omit()
non.dune_ilr_df <- subset(diff_iso_ilr_df, ecotype=="non.dune") %>% 
  na.omit()

# split data frames into list of dataframes so we can do the ilr transform
# for each gene
dune_ilr_list <- split( dune_ilr_df , f = dune_ilr_df$gene )
non.dune_ilr_list <- split( non.dune_ilr_df, f = non.dune_ilr_df$gene)

dune_res <- list()
i <- 1
for( gene_df in dune_ilr_list ){
  info <- gene_df[1:3]
  info$ilr <- as.vector(ilr(gene_df[,4:5]))
  dune_res[[i]] <- info
  i <- i + 1
}

non.dune_res <- list()
i <- 1
for( gene_df in  non.dune_ilr_list ){
  info <- gene_df[1:3]
  info$ilr <- as.vector(ilr(gene_df[,4:5]))
  non.dune_res[[i]] <- info
  i <- i + 1
}

diff_iso_splice_df <- bind_rows(dune_res, non.dune_res) %>%
  dplyr::select(!ecotype) %>%
  pivot_wider(names_from=sample, values_from = ilr) %>%
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
  geom_point(size=4, alpha=.75) +
  geom_polygon(data=hull_diff_iso_d, alpha=0.5) +
  geom_polygon(data=hull_diff_iso_nd, alpha=0.5) +
  
  labs(x=paste0("PC1 ","(",round(diff_iso_splice_pve_df[1,2], 2),"%)"),
       y=paste0("PC2 ", "(",round(diff_iso_splice_pve_df[2,2],2),"%)")) +#,
  #title="Alternative splicing") +
  theme_bw(base_size=12) +
  theme(legend.position = "none",
        legend.title = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  geom_hline(yintercept = 0, lty = 2, col = "grey50") +
  geom_vline(xintercept = 0, lty = 2, col = "grey50") +
  scale_color_manual(values=c("gold2", "forestgreen")) +
  scale_fill_manual(values=c("gold2", "forestgreen")) +
  scale_shape_manual(values=c(24,21))

ggsave("figures/FIG_Sx_diff_iso_PCA.pdf", plot=plot_diff_iso_splice_pca,
       device = "pdf", height = 65.625 , width = 87.5, dpi = 300, units = "mm")

#### plot splicing of individual genes ####
library(ggridges)

plot_GLH17_iso_props <- ggplot(diff_iso_all_prop_df %>%
                                 filter(gene=="TRINITY_DN6261_c0_g1.0"),
                               aes(x = proportion_TPM, y = fct_rev(as.factor(isoform_id)), fill=ecotype)) +
  geom_density_ridges(scale=1, alpha = .5) +
  facet_wrap(~ecotype) +
  theme_bw(base_size=18) +
  theme(legend.position = "none") +
  scale_fill_manual(values=c("gold2", "forestgreen")) +
  labs(y="", x="Isoform proportion")

ggsave("figures/plot_GLH17_raw.png", plot = (plot_spacer() | plot_GLH17_iso_props) / (plot_GLH17_expr | plot_GLH17_psi),
       device = "png", height = 9, width = 12, dpi = 300, units = "in")

plot_iso_props <-  ggplot(diff_iso_all_prop_df %>%
                                    filter(gene=="TRINITY_DN6473_c1_g1.0"),
                                  aes(x = proportion_TPM, y = fct_rev(as.factor(isoform_id)), fill=ecotype)) +
  geom_density_ridges(scale=1, alpha = .5) +
  facet_wrap(~ecotype) +
  theme_bw(base_size=18) +
  theme(legend.position = "none",
        axis.text.y = element_blank()) +
  scale_fill_manual(values=c("gold2", "forestgreen")) +
  labs(y="", x="Isoform proportion")


plot_ABC1K14_iso_props <-  ggplot(diff_iso_all_prop_df %>%
                                    filter(gene=="TRINITY_DN45768_c0_g1.0"),
                                  aes(x = proportion_TPM, y = fct_rev(as.factor(isoform_id)), fill=ecotype)) +
  geom_density_ridges(scale=1, alpha = .5) +
  facet_wrap(~ecotype) +
  theme_bw(base_size=18) +
  theme(legend.position = "none") +
  scale_fill_manual(values=c("gold2", "forestgreen")) +
  labs(y="", x="Isoform proportion")

#TRINITY_DN51265_c0_g1
plot_LHCB5_iso_props <- ggplot(diff_iso_all_prop_df %>%
                                 filter(gene=="TRINITY_DN51265_c0_g1.0"),
                               aes(x = proportion_TPM, y = fct_rev(as.factor(isoform_id)), fill=ecotype)) +
  geom_density_ridges(scale=1, alpha = .5) +
  facet_wrap(~ecotype) +
  theme_bw(base_size=18) +
  theme(legend.position = "none") +
  scale_fill_manual(values=c("gold2", "forestgreen")) +
  labs(y="", x="Isoform proportion")

# rhodanese TRINITY_DN6462_c0_g1 Ha412HOChr10g0445261
plot_DN6462_c0_g1_iso_props <- ggplot(diff_iso_all_prop_df %>%
                  filter(gene=="TRINITY_DN6462_c0_g1.0"),
                aes(x = proportion_TPM, y = fct_rev(as.factor(isoform_id)), fill=ecotype)) +
  geom_density_ridges(scale=1, alpha = .5) +
  facet_wrap(~ecotype) +
  theme_bw(base_size=18) +
  theme(legend.position = "none") +
  scale_fill_manual(values=c("gold2", "forestgreen")) +
  labs(y="", x="Isoform proportion")

# AT1G08520. Encodes the CHLD subunit of the Mg-chelatase enzyme involved in chlorophyll biosynthesis. Lines carrying recessive mutations of this locus are white and seedling lethal.
plot_ALB1_iso_props <- ggplot(diff_iso_all_prop_df %>%
                                 filter(gene=="TRINITY_DN20639_c0_g1.0"),
                               aes(x = proportion_TPM, y = fct_rev(as.factor(isoform_id)), fill=ecotype)) +
  geom_density_ridges(scale=1, alpha = .5) +
  facet_wrap(~ecotype) +
  theme_bw(base_size=18) +
  theme(legend.position = "none") +
  scale_fill_manual(values=c("gold2", "forestgreen")) +
  labs(y="", x="Isoform proportion")

# Ha412HOChr11g0508211. AT2G42450
# COR15A/LEA24? I think mislabeled on Arabidopsis.org? not going to highlight this one
plot_COR15A_iso_props <- ggplot(diff_iso_all_prop_df %>%
                               filter(gene=="TRINITY_DN11567_c1_g1.0"),
                               aes(x = proportion_TPM, y = fct_rev(as.factor(isoform_id)), fill=ecotype)) +
  geom_density_ridges(scale=1, alpha = .5) +
  facet_wrap(~ecotype) +
  theme_bw(base_size=18) +
  theme(legend.position = "none") +
  scale_fill_manual(values=c("gold2", "forestgreen")) +
  labs(y="", x="Isoform proportion")

plot_SCPL13_iso_props <- ggplot(diff_iso_all_prop_df %>%
                                  filter(gene=="TRINITY_DN11324_c1_g1.0"),
                                aes(x = proportion_TPM, y = isoform_id, fill=ecotype)) +
  geom_density_ridges(scale = 1, alpha = .5) +
  facet_wrap(~ecotype) +
  theme_bw(base_size=18) +
  theme(legend.position = "none") +
  scale_fill_manual(values=c("gold2", "forestgreen")) +
  labs(y="", x="Isoform proportion")

plot_CESA6_iso_props <- ggplot(diff_iso_all_prop_df %>%
                                 filter(gene=="TRINITY_DN3311_c0_g1.0"),
                               aes(x = proportion_TPM, y = isoform_id, fill=ecotype)) +
  geom_density_ridges(scale=1, alpha = .5) +
  facet_wrap(~ecotype) +
  theme_bw(base_size=18) +
  theme(legend.position = "none") +
  scale_fill_manual(values=c("gold2", "forestgreen")) +
  labs(y="", x="Isoform proportion")

#TRINITY_DN22648_c0_g1; EMB2004; AT1G10510
plot_EMB2004_iso_props <- ggplot(diff_iso_all_prop_df %>%
                                 filter(gene=="TRINITY_DN22648_c0_g1.0"),
                               aes(x = proportion_TPM, y = isoform_id, fill=ecotype)) +
  geom_density_ridges(scale=1, alpha = .5) +
  facet_wrap(~ecotype) +
  theme_bw(base_size=18) +
  theme(legend.position = "none") +
  scale_fill_manual(values=c("gold2", "forestgreen")) +
  labs(y="", x="Isoform proportion")

#### trying just PSI instead of ILR ####
tmp <- dplyr::select(diff_iso_ilr_df, gene, sample, iso_1) %>% pivot_wider(names_from = sample, values_from = iso_1) %>%
  dplyr::select(!gene) %>%
  na.omit()

diff_iso_splice_pca <- rda(t(tmp), scale = F, center=T)
            