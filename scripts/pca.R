#### SNP PCA ####

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


# convert VCF file to GDS format
vcf.fn <- "data/hc_out/dune_non-dune.nHet_filtered.hc.vcf.gz"
snpgdsVCF2GDS(vcf.fn, "analysis/pca/vcf.gds", method = "biallelic.only")
snpgdsSummary("analysis/pca/vcf.gds")

#### Analysis ####

# population information from text file
pop_code <- scan("analysis/pca/pops.txt", what=character()) # "pops.txt" should be a single line with "population" noted for each sample, in quotes, tab separated. order needs to match order of samples in the vcf.

# open the GDS file
genofile <- snpgdsOpen("analysis/pca/vcf.gds")

# Prune for LD. 
?snpgdsLDpruning
set.seed(1000)
snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2, slide.max.bp = 50000, autosome.only = FALSE) # Try different LD thresholds for sensitivity analysis? this step takes a while
save(snpset, file="analysis/pca/snpset.list")
load("~/gsd_RNA-seq/analysis/pca/snpset.list")

# Get all selected snp id
snpset.id <- unlist(snpset)

# Run PCA. 
snp_pca <- snpgdsPCA(genofile, num.thread=2, autosome.only = FALSE)
snp_pca_LD_pruned <- snpgdsPCA(genofile, snp.id=snpset.id, num.thread=2, autosome.only = FALSE)
head(snp_pca)

# Get sample ids, make sure they match up with population codes
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
#View((cbind(sample.id, pop_code))) #check that they match

# close the genotype file
snpgdsClose(genofile)

# Make a data.frame
snp_pca_df <- data.frame(sample.id = snp_pca$sample.id,
                  population = factor(pop_code)[match(snp_pca$sample.id, sample.id)],
                  EV1 = snp_pca$eigenvect[,1],    # the first eigenvector
                  EV2 = snp_pca$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)
#head(snp_pca_df)


# variance proportion (% variance explained)
snp_pve_df <- data.frame(pc = 1:24, pve = snp_pca$varprop*100)

plot_snp_pve <- ggplot(snp_pve_df, aes(x=pc, y=pve)) +
  geom_bar(stat = "identity") +
  xlab("PC") +
  ylab("Percent variance explained")
  
# plot in ggplot
plot_snp_pca <- ggplot(snp_pca_df, aes(x=EV1, y=EV2, color=population, shape=population)) + 
  geom_point(size=3, alpha=.5) +
  labs(x=paste0("PC1 ","(",round(snp_pve_df[1,2], 2),"%)"), y=paste0("PC2 ", "(",round(snp_pve_df[2,2],2),"%)"), title="SNPs") +
  theme_bw() +
  theme(text = element_text(size=14),
        legend.title=element_blank(),
        legend.position = "bottom",
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,-10,-10,-10)) +
  scale_color_manual(values=c("gold2", "forestgreen")) +
  scale_shape_manual(values=c(17,19))
  

#### EXPRESSION PCA ####

#BiocManager::install('PCAtools')

# PCA with 'regularized log'-transformed gene expression levels from DESeq2 (see analyze_differential_expression_DESeq2.R script for details on how we get the object rlog_counts.)
sample_table <- read.table("analysis/DESeq2/sample_table.txt", header=T)
expr_pca <- rda(t(rlog_counts), scale = F, center=T) #have to transpose the dataframe
str(expr_pca)
summary(expr_pca)

# variance proportion (% variance explained)
expr_pve_df <- data.frame(summary(eigenvals(expr_pca)))[2,1:12] %>%
  pivot_longer(values_to = "PVE", names_to="PC", cols = 1:12) %>%
  mutate(PVE=100*PVE)

#expr_pca_df <- data.frame(expr_pca_prcomp$x) # pull out the PCs; they are in the matrix 'x' within the expr_pca object
#expr_pca_df <- tibble::rownames_to_column(expr_pca_df, "sample.id") # convert rownames (sample_ids) to the first column

#expr_pop_codes <- read.table("expr_pops.txt", header = T)
#expr_pop_codes <- read.table("expr_pops.ebseq.txt", header=T)

#expr_pca_df$population <- as.factor(expr_pop_codes$population)
#expr_pca_df <- select(expr_pca_df, sample.id, population, PC1, PC2, PC3)
#str(expr_pca_df)

#### plotting prcomp PCA ####
#plot_expr_pca <- ggplot(expr_pca_df, aes(x=PC1, y=PC2, color=population)) + 
#  geom_point(size=3, alpha=.5) +
#  labs(x=paste0("PC1 ","(",round(expr_pve_df[1,2], 2),"%)"), y=paste0("PC2 ", "(",round(expr_pve_df[2,2],2),"%)"),
#       title="Gene expression") +
#  theme_bw() +
#  theme(text = element_text(size=14),
#        legend.title=element_blank(),
#        legend.position = "none") +
#  scale_color_manual(values=c("gold2", "forestgreen"))

# get site scores for plotting
populations <- c("d2", "d3", "d3", "d3", "d3", "nd1", "nd1", "nd1", "nd1", "d1", "nd2", "nd2", "nd2", "nd2", "nd3", "nd3", "nd3", "d1", "nd3", "d1", "d1", "d2", "d2", "d2") # to check if samples cluster by sampling location within habitat (they don't)
expr_pca.site_sc <- data.frame(scores(expr_pca, choices = 1:2, scaling = 1, display = "wa")) %>%
  tibble::rownames_to_column("sample") %>%
  mutate(habitat=sample_table$habitat) %>%
  mutate(population=populations)

#expr_pca_prcomp.site_sc <- data.frame(scores(expr_pca_prcomp, choices = 1:2, display = "sites"))

# plotting
plot_expr_pca <-  ggplot(data=expr_pca.site_sc, aes(x = PC1, y = PC2, color=habitat, shape=habitat)) +
  geom_point(size=3, alpha=.5) +
  labs(x=paste0("PC1 ","(",round(expr_pve_df[1,2], 2),"%)"), y=paste0("PC2 ", "(",round(expr_pve_df[2,2],2),"%)"),
       title="Gene expression") +
  theme_bw() +
  theme(text = element_text(size=14),
        legend.position = "none",
        legend.title = element_blank()) +
  scale_color_manual(values=c("gold2", "forestgreen")) +
  scale_shape_manual(values=c(17,19))

#### splicing PCA? ####
# Using the "% Spliced In" values from rMATS output. In this data set, each row is an AS event, each column is a plant, so need to transpose before PCA
splice_pca_df <- all_AS_events_PSI[,4:27] %>%
  lapply(as.numeric) %>% as.data.frame() %>%
  na.omit() 
splice_pca <- rda(t(splice_pca_df), scale = F, center=T)

splice_pve_df <- data.frame(summary(eigenvals(splice_pca)))[2,1:12] %>%
  pivot_longer(values_to = "PVE", names_to="PC", cols = 1:12) %>%
  mutate(PVE=100*PVE)

sample_list <- read.table("sample_list.txt", header=T)
splice_pca.site_sc <- data.frame(scores(splice_pca, choices = 1:2, scaling = 1, display = "wa")) %>%
  tibble::rownames_to_column("sample") %>%
  mutate(habitat=sample_list$habitat)

# plotting
plot_splice_pca <-  ggplot(data=splice_pca.site_sc, aes(x = PC1, y = PC2, color=habitat, shape=habitat)) +
  geom_point(size=3, alpha=.5) +
  labs(x=paste0("PC1 ","(",round(splice_pve_df[1,2], 2),"%)"), y=paste0("PC2 ", "(",round(splice_pve_df[2,2],2),"%)"),
       title="Splicing") +
  theme_bw() +
  theme(text = element_text(size=14),
        legend.position = "none",
        legend.title = element_blank()) +
  scale_color_manual(values=c("gold2", "forestgreen")) +
  scale_shape_manual(values=c(17,19))


#### join plots ####
#FIG_1 <- plot_grid(plot_snp_pca, plot_expr_pca, labels=c("a)","b)"), ncol=1, nrow=2) #using cowplot package
patchwork <- gsd_sampling_map | (plot_snp_pca / plot_expr_pca) 
FIG_1 <- patchwork + plot_annotation(tag_levels = "A") #gsd_sampling_map is from sampling_map.R script

#arrange with patchwork
ggsave("figures/FIG_1.jpg", plot=FIG_1, width=17, dpi = 600, units = "cm")

#### k-means clustering?? ####
?kmeans
