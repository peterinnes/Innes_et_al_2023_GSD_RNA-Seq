# Perform PCA on SNP, expression, and splicing data

#if (!requireNamespace("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")
#install.packages("devtools")
#BiocManager::install("SNPRelate")
#devtools::install_github("thomasp85/patchwork")

library(SNPRelate)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(vegan)

#### SNP pca ####

# convert VCF file to GDS format
showfile.gds(closeall=TRUE)
vcf.fn <- "data/hc_out/dune_non-dune.nHet_filtered.hc.vcf.gz"
vcf_pruned.fn <- "data/hc_out/plink/pruneddata.vcf" #pruned with PLINK, LD < .2, 500kb windows
snpgdsVCF2GDS(vcf.fn, "analysis/pca/vcf.gds")
snpgdsVCF2GDS(vcf_pruned.fn, "analysis/pca/vcf_pruned.gds")
#snpgdsSummary("analysis/pca/vcf.gds")
snpgdsSummary("analysis/pca/vcf_pruned.gds")

# habitat information from text file
pop_code <- scan("analysis/pca/pops.txt", what=character()) # "pops.txt" should be a single line with "habitat" i.e. ecotype noted for each sample, in quotes, tab separated. order needs to match order of samples in the vcf.
habitat_pop_sample <- read.table("analysis/pca/habitat_pop_sample.txt", header=T) #this file also has the specific sampling site aka population, 'pop' for each sample

# open the GDS file
genofile <- snpgdsOpen("analysis/pca/vcf.gds")
genofile <- snpgdsOpen("analysis/pca/vcf_pruned.gds")

# Prune for LD with SNPrelate. This uses slightly different algo compared to PLINK, I think
?snpgdsLDpruning
set.seed(1000)
snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2,
                          method = "r", autosome.only = F)
##Try different LD thresholds for sensitivity analysis? this step takes a while
#save(snpset, file="analysis/pca/snpset.list")
#load("~/gsd_RNA-seq/analysis/pca/snpset.list")
# Get all selected snp id
snpset.id <- unlist(snpset)

# Run PCA. 
snp_pca <- snpgdsPCA(genofile, num.thread=2, autosome.only = FALSE)
snp_pca <- snpgdsPCA(genofile, snp.id=snpset.id,
                     num.thread=2, autosome.only = FALSE) #only for SNPrelate LD-pruned
head(snp_pca)

# Get sample ids, make sure they match up with habitat codes
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
#View((cbind(sample.id, pop_code))) #check that they match

# close the genotype file
snpgdsClose(genofile)

# Make a data.frame
snp_pca_df <- data.frame(sample.id = snp_pca$sample.id,
                  habitat = factor(pop_code)[match(snp_pca$sample.id, sample.id)],
                  EV1 = snp_pca$eigenvect[,1],    # the first eigenvector
                  EV2 = snp_pca$eigenvect[,2],    # the second eigenvector
                  EV3 = snp_pca$eigenvect[,3],
                  EV4 = snp_pca$eigenvect[,4],
                  EV5 = snp_pca$eigenvect[,5],
                  stringsAsFactors = FALSE)

snp_pca_df <- left_join(snp_pca_df, habitat_pop_sample) %>%
  rename(habitat="ecotype")
#head(snp_pca_df)


# variance proportion (% variance explained)
snp_pve_df <- data.frame(pc = 1:24, pve = snp_pca$varprop*100)

plot_snp_pve <- ggplot(snp_pve_df, aes(x=pc, y=pve)) +
  geom_bar(stat = "identity") +
  xlab("PC") +
  ylab("Percent variance explained")

# get convex hulls for pca plot
hull_snp_d <- filter(snp_pca_df, ecotype=="dune") %>%
  dplyr::slice(chull(EV1,EV2))
hull_snp_nd <- filter(snp_pca_df, ecotype=="non-dune") %>%
  dplyr::slice(chull(EV1,EV2))

# plot in ggplot
plot_snp_pca <- ggplot(snp_pca_df,
                                 aes(x=EV1, y=EV2, shape=ecotype, fill=ecotype)) +
  # plot each plant and plot convex hulls for each ecotype
  geom_point(size=5, alpha=.75) +
  geom_polygon(data=hull_snp_d, alpha=.5) +
  geom_polygon(data=hull_snp_nd, alpha=.5) +
  
  labs(x=paste0("PC1 ","(",round(snp_pve_df[1,2], 2),"%)"),
       y=paste0("PC2 ", "(",round(snp_pve_df[2,2],2),"%)")) + #, title="SNPs") +
  theme_bw() +
  theme(text = element_text(size=24),#,
        #legend.title=element_blank(),
        legend.position = "none",
        #legend.position = c(.75,.75),
        #legend.background = element_blank(),
        #legend.box.background = element_rect(color = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.title = element_text(size=18),
        axis.text = element_text(size=18)) +
  geom_hline(yintercept = 0, lty = 2, col = "grey50") +
  geom_vline(xintercept = 0, lty = 2, col = "grey50") +
        #legend.margin=margin(0,0,0,0),
        #legend.box.margin=margin(0,0,0,0)) +
  #scale_color_manual(values=c("gold2", "forestgreen")) +
  scale_fill_manual(values = c("gold2", "forestgreen")) +
  scale_shape_manual(values=c(24,21))
plot_snp_pca

plot_snp_pca_5v4 <- ggplot(snp_pca_df,
                                     aes(x=EV4, y=EV5,
                                         color=ecotype, shape=ecotype)) + 
  geom_point(size=5, alpha=.5) +
  labs(x=paste0("PC4 ","(",round(snp_pve_df[4,2], 2),"%)"),
       y=paste0("PC5 ", "(",round(snp_pve_df[5,2],2),"%)"), title="SNPs") +
  theme_bw() +
  theme(text = element_text(size=18),#,
        #legend.title=element_blank(),
        legend.position = "none") +
  #legend.margin=margin(0,0,0,0),
  #legend.box.margin=margin(0,0,0,0)) +
  scale_color_manual(values=c("gold2", "forestgreen")) +
  scale_shape_manual(values=c(24,21))

#### EXPRESSION pca ####

#BiocManager::install('PCAtools')

# PCA with 'regularized log'-transformed gene expression levels from DESeq2
# (see analyze_differential_expression_DESeq2.R script for details 
# on how we get the object rlog_counts.)
sample_table_deseq2 <- read.table("analysis/DESeq2/sample_table.txt", header=T)
rlog_counts <- read.csv("analysis/pca/rlog-transformed_expression_counts.csv")

#have to transpose the dataframe. scale=F because rlog_counts already on similar scale.
expr_pca <- rda(t(rlog_counts), scale = F, center=T) 
summary(expr_pca)

# percent variance explained
expr_pve_df <- data.frame(summary(eigenvals(expr_pca)))[2,1:12] %>%
  pivot_longer(values_to = "PVE", names_to="PC", cols = 1:12) %>%
  mutate(PVE=100*PVE)

# get PC site scores for plotting
expr_pca.site_sc <- data.frame(scores(expr_pca,
                                      choices = 1:2, scaling = 1, display = "wa")) %>%
  tibble::rownames_to_column("sample.id") %>%
  mutate(habitat=sample_table_deseq2$habitat)# %>%
  #mutate(population=populations)

# get convex hulls for expression pca plot
hull_expr_d <- filter(expr_pca.site_sc, habitat=="dune") %>%
  dplyr::slice(chull(PC1,PC2))
hull_expr_nd <- filter(expr_pca.site_sc, habitat=="non-dune") %>%
  dplyr::slice(chull(PC1,PC2))

# plotting
plot_expr_pca <-  ggplot(data=expr_pca.site_sc, aes(x = PC1, y = PC2,
                                                    shape=habitat,
                                                    fill=habitat)) +
  geom_point(size=5, alpha=.75) +
  geom_polygon(data=hull_expr_d, alpha=.5) +
  geom_polygon(data=hull_expr_nd, alpha=.5) +
  
  labs(x=paste0("PC1 ","(",round(expr_pve_df[1,2], 2),"%)"),
       y=paste0("PC2 ", "(",round(expr_pve_df[2,2],2),"%)")) + #,
       #title="Transcript levels") +
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
        #legend.position = "bottom",
        #legend.margin=margin(0,0,0,0),
        #legend.box.margin=margin(0,0,0,0)) +
  #scale_color_manual(values=c("gold2", "forestgreen")) +
  scale_fill_manual(values=c("gold2", "forestgreen")) +
  scale_shape_manual(values=c(24,21))

#### SPLICING pca ####
# Using the "% Spliced In" i.e. PSI values from rMATS output. In this data set, each row is an AS event, each column is a plant, so need to transpose before PCA. Getting rid of rows with any missing values takes us down to 16065 splice events, from 18038
sample_list <- read.table("sample_list.txt", header=T)

rmats_splice_pca_df <- all_AS_events_PSI[,4:27] %>%
  lapply(as.numeric) %>% as.data.frame() #%>% na.omit()
rownames(rmats_splice_pca_df) <- all_AS_events_PSI$ID

# remove events with > 40% missingness
rmats_splice_pca_df <- rmats_splice_pca_df %>%
  mutate(missing_perc = rowMeans(is.na(rmats_splice_pca_df))) %>%
  filter(missing_perc<.4) %>%
  dplyr::select(!missing_perc)

# set missing values to the average PSI of each event (each row)
for( i in 1:nrow(rmats_splice_pca_df) ){
  for( j in 1:ncol(rmats_splice_pca_df) ){
    if( is.na(rmats_splice_pca_df[i,j]) ){
      rmats_splice_pca_df[i,j] <- rowMeans(rmats_splice_pca_df[i,], na.rm = T)
    }
  }
}

rmats_splice_pca <- rda(t(rmats_splice_pca_df), scale = F, center=T)

## compare to splice PSI values that have undergone standardization
## and quantile normalization
#rmats_splice_pca <- rda(t(rmats_splice_pheno_df), scale=F, center=F)

rmats_splice_pve_df <- data.frame(summary(
  eigenvals(rmats_splice_pca)))[2,1:12] %>%
  pivot_longer(values_to = "PVE",
               names_to="PC",
               cols = 1:12) %>%
  mutate(PVE=100*PVE)

rmats_splice_pca.site_sc <- data.frame(scores(rmats_splice_pca, choices = 1:2,
                                              scaling = 1, display = "wa")) %>%
  tibble::rownames_to_column("sample.id") %>%
  mutate(habitat=sample_list$habitat) %>%
  mutate(sample.id=gsub("_R1_001.trimmed.fq.gz","",sample.id),
         sample.id=gsub("\\.","-",sample.id))

# get convex hulls for splicing pca plot
hull_rmats_d <- filter(rmats_splice_pca.site_sc, habitat=="dune") %>%
  dplyr::slice(chull(PC1,PC2))
hull_rmats_nd <- filter(rmats_splice_pca.site_sc, habitat=="non-dune") %>%
  dplyr::slice(chull(PC1,PC2))

# plotting
plot_rmats_splice_pca <-  ggplot(data=rmats_splice_pca.site_sc,
                                 aes(x = PC1, y = PC2, shape=habitat,
                                     fill=habitat)) +
  geom_point(size=5, alpha=.75) +
  geom_polygon(data=hull_rmats_d, alpha=0.5) +
  geom_polygon(data=hull_rmats_nd, alpha=0.5) +
  
  labs(x=paste0("PC1 ","(",round(rmats_splice_pve_df[1,2], 2),"%)"),
       y=paste0("PC2 ", "(",round(rmats_splice_pve_df[2,2],2),"%)")) +#,
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

# DEXSeq splicing (DEU) pca
sample_table_dexseq <- read.table("sample_list.txt", header=T)
rlog_exon_counts <- read.csv("analysis/pca/dexseq_exon_counts.csv")
dexseq_splice_pca <-rda(t(rlog_exon_counts), scale=F, center=T) 

dexseq_splice_pve_df <- data.frame(summary(
  eigenvals(dexseq_splice_pca)))[2,1:12] %>%
  pivot_longer(values_to = "PVE",
               names_to="PC",
               cols = 1:12) %>%
  mutate(PVE=100*PVE)

dexseq_splice_pca.site_sc <- data.frame(scores(dexseq_splice_pca,
                                               choices = 1:4,
                                               scaling = 1,
                                               display = "wa")) %>%
  tibble::rownames_to_column("sample.id") %>%
  mutate(habitat=sample_table_dexseq$habitat)

plot_dexseq_splice_pca <-  ggplot(data=dexseq_splice_pca.site_sc,
                                  aes(x = PC1, y = PC2,
                                      color=habitat, shape=habitat)) +
  geom_point(size=5, alpha=.5) +
  labs(x=paste0("PC1 ","(",round(dexseq_splice_pve_df[1,2], 2),"%)"),
       y=paste0("PC2 ", "(",round(dexseq_splice_pve_df[2,2],2),"%)"),
       title="Splicing (exon-based)") +
  theme_bw() +
  theme(text = element_text(size=18),
        #legend.position = "none",
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(0,0,0,0),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.title = element_text(size=18),
        axis.text = element_text(size=12)) +
  geom_hline(yintercept = 0, lty = 2, col = "grey50") +
  geom_vline(xintercept = 0, lty = 2, col = "grey50") +
  scale_color_manual(values=c("gold2", "forestgreen")) +
  scale_shape_manual(values=c(24,21))

# leafcutter PCA (intron excision ratios)
sample_table_leafcutter <- read.table("data2/leafcutter_out/groups_file.txt", col.names=c("sample.id", "habitat")) %>%
  mutate(sample.id=gsub("-",".", sample.id))

leafcutter_splice_pca <- rda( t(ier_df[-c(1:4)]),scale=F, center=F)

leafcutter_splice_pve_df <- data.frame(summary(eigenvals(leafcutter_splice_pca)))[2,1:12] %>%
  pivot_longer(values_to = "PVE", names_to="PC", cols = 1:12) %>%
  mutate(PVE=100*PVE)

leafcutter_splice_pca.site_sc <- data.frame(scores(leafcutter_splice_pca, choices = 1:4, scaling = 1, display = "wa")) %>%
  tibble::rownames_to_column("sample.id") %>%
  full_join(sample_table_leafcutter)

plot_leafcutter_splice_pca <-  ggplot(data=leafcutter_splice_pca.site_sc, aes(x = PC1, y = PC2, color=habitat, shape=habitat)) +
  geom_point(size=5, alpha=.5) +
  labs(x=paste0("PC1 ","(",round(leafcutter_splice_pve_df[1,2], 2),"%)"), y=paste0("PC2 ", "(",round(leafcutter_splice_pve_df[2,2],2),"%)"),
       title="Splicing (intron excision)") +
  theme_bw() +
  theme(text = element_text(size=18),
        #legend.position = "none",
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(0,0,0,0),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.title = element_text(size=18),
        axis.text = element_text(size=12)) +
  geom_hline(yintercept = 0, lty = 2, col = "grey50") +
  geom_vline(xintercept = 0, lty = 2, col = "grey50") +
  scale_color_manual(values=c("gold2", "forestgreen")) +
  scale_shape_manual(values=c(24,21))


#### join plots ####
#FIG_1 <- plot_grid(plot_snp_pca, plot_expr_pca, labels=c("a)","b)"), ncol=1, nrow=2) #using cowplot package
FIG_1 <- (plot_spacer() | plot_snp_pca) / (plot_expr_pca | plot_rmats_splice_pca)
#FIG_1 <- patchwork + plot_annotation(tag_levels = "b") #gsd_sampling_map is from sampling_map.R script

#poster_patchwork <- (plot_snp_pca / plot_expr_pca) | (plot_dexseq_splice_pca / plot_rmats_splice_pca)
#poster_FIG_1 <- poster_patchwork + plot_annotation(tag_levels = "A")

#arrange with patchwork
ggsave("figures/FIG_1_raw.png", plot=FIG_1,device = "png",
       height = 9, width = 12, dpi = 300, units = "in")

#ggsave("figures/poster_FIG_1", plot=poster_patchwork, device = "png", dpi=300, width=8, height=6, units="in")

ggsave("figures/SNP_PCA.png", plot=plot_snp_pca, device = "png", dpi=300, width=8, height=6, units="in")
ggsave("figures/expr_PCA.png", plot=plot_expr_pca, device = "png", dpi=300, width=8, height=6, units="in")
ggsave("figures/DEU_PCA.png", plot=plot_dexseq_splice_pca, device = "png", dpi=300, width=8, height=6, units="in")
ggsave("figures/rMATS_PCA.png", plot=plot_rmats_splice_pca, device = "png", dpi=300, width=8, height=6, units="in")
#### k-means clustering?? ####
?kmeans
