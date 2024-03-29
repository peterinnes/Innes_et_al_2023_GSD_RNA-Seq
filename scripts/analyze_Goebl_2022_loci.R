library(GenomicRanges)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(reshape2)

# get ranges of all expressed genes
expressed_genes <- read.table("data/expressed_genes.txt",
                              col.names = "Ha412_gene")
genes_gff <- read.table("data/HAN412_Eugene_curated_v1_1.genes_gff.tmp") %>% 
  dplyr::select(c(1:5,11)) %>%
  filter(!grepl("Chr00", V1)) #exclude ~550 genes on unplaced contigs

names(genes_gff) <- c("chrom","gene_start", "gene_end",
                      "gene_length", "strand", "Ha412_gene")
expressed_genes_gff <- genes_gff %>% filter(Ha412_gene %in% expressed_genes$Ha412_gene)
expressed_genes_ranges <- makeGRangesFromDataFrame(expressed_genes_gff,
                                                   keep.extra.columns = T)

# DS genes
DS_genes <- read.table("data/study_DS_rMATS_genes.txt",
                       col.names = c("Ha412_gene"))

# DE + DS genes
DE_DS_union_genes <- read.table("data/study_DE_genes_noLFCthreshold.txt") %>%
  dplyr::rename(Ha412_gene=V1) %>%
  full_join(DS_genes) %>%
  unique()

# non-DE/DS genes i.e. "control" genes
control_genes <- expressed_genes_gff %>%
  filter(!Ha412_gene %in% DE_DS_union_genes$Ha412_gene) %>%
  dplyr::select(Ha412_gene)

# get loci from Goebl et al 2022 
goebl_2022_loci_95_quant <- read.csv("data/20230609_GSDseqAnalyCombDat.csv",
                header = T) %>% 
    dplyr::rename(chrom=CHROM, pos=POS) %>%
    filter(DIFFS_Hd_ABS >= quantile(DIFFS_Hd_ABS, .95)) 

goebl_2022_loci_99_quant <- read.csv("data/20230609_GSDseqAnalyCombDat.csv",
                                  header = T) %>%
    dplyr::rename(chrom=CHROM, pos=POS) %>%
    filter(DIFFS_Hd_ABS >= quantile(DIFFS_Hd_ABS, .99))

goebl_2022_ranges_95_quant <- makeGRangesFromDataFrame(goebl_2022_loci_95_quant %>%
                                                           mutate(start=pos-1, end=pos) %>%
                                                           dplyr::select(chrom, start, end, DIFFS_Hd_ABS),
                                                       keep.extra.columns = T)

goebl_2022_ranges_99_quant <- makeGRangesFromDataFrame(goebl_2022_loci_99_quant %>%
                                                           mutate(start=pos-1, end=pos) %>%
                                                           dplyr::select(chrom, start, end, DIFFS_Hd_ABS),
                                                       keep.extra.columns = T)

#### get distance between genes of interest and Goebl et al 2022 SNPs ####
# code below is adapted from Verta & Jones 2019 eLife. 

# get distance between all expressed genes and closest locus adaptive locus
# IMPORTANT: have to update 'hits' object 
# for 95 quantile threshold vs 99 quantile threshold and then re-run subsequent code.
# didn't have time to make a loop...
expressed_genes_ranges$dist_to_adaptive <- NA
hits <- distanceToNearest(expressed_genes_ranges, goebl_2022_ranges_95_quant)
expressed_genes_ranges$dist_to_adaptive[queryHits(hits)] <- mcols(hits)$distance

# split into test and control groups
test_ranges <- expressed_genes_ranges[which(expressed_genes_ranges$Ha412_gene %in%
                                          DS_genes$Ha412_gene)]
control_ranges <- expressed_genes_ranges[which(expressed_genes_ranges$Ha412_gene %in%
                                             control_genes$Ha412_gene)]

# This loop calculates the proportion of DS genes and control genes 
# x bp or closer to a Goebl et al 2022 SNP,  where x is the current distance ("window"). 
# Each iteration increases the distance from the gene.
windows <- seq(1,5000000,50000)
dist_to_adaptive <- data.frame()
for (i in 1:length(windows)){
  dist_to_adaptive[i,'test'] <- length(test_ranges[which(test_ranges$dist_to_adaptive < windows[i])]$Ha412_gene) /
    length(test_ranges$Ha412_gene)
  dist_to_adaptive[i,'control'] <- length(control_ranges[which(control_ranges$dist_to_adaptive < windows[i])]$Ha412_gene) /
    length(control_ranges$Ha412_gene)
  dist_to_adaptive[i,'distance'] <- windows[i]
}
#rownames(dist_to_adaptive) = as.character(windows)
# reshape the df for plotting
dist_to_adaptive <- dist_to_adaptive %>%
  pivot_longer(1:2, names_to = "set", values_to = "proportion")
dist_to_adaptive$set <- factor(dist_to_adaptive$set, levels = c("test", "control"))

#### random expectation ####
random_expectation_dist = data.frame()
for (y in 1:100){
  for (i in 1:length(windows)){
    r <- data.frame(control_ranges) %>%
      filter(Ha412_gene %in% sample(control_ranges$Ha412_gene,
                                    length(test_ranges$Ha412_gene)))
    
    random_expectation_dist[i,y] <- nrow(
      subset(r, dist_to_adaptive < windows[i])) / length(r$Ha412_gene)
  }
}

rownames(random_expectation_dist) = as.character(windows)

quantile_95 = function(x){quantile(x,probs=c(0.05,0.95))}
random_expectation_95CI = apply(random_expectation_dist, 1, quantile_95)
random_expectation_95CI = data.frame(t(random_expectation_95CI)) %>%
  dplyr::rename(min=X5., max=X95.) %>%
  rownames_to_column("distance")
random_expectation_95CI$distance <- as.numeric(random_expectation_95CI$distance)

#### plot ####
plot_dist_to_goebl_loci_95_quant <- ggplot(data = dist_to_adaptive) +
  geom_line(aes(x=distance, y=proportion, linetype=set, color=set),
            linewidth=.25) +
  scale_color_manual(values = c("red", "black")) +
  scale_linetype_manual(values = c(1,2)) +
  geom_ribbon(data=random_expectation_95CI,
              aes(x=distance, ymin = min, ymax = max),
              color='grey50', alpha=.25, linewidth=.25) +
  scale_x_continuous(limits = c(0,1000001)) +
  scale_y_continuous(limits = c(0,.3)) +
  theme_bw(base_size = 12) +
  theme(legend.position =  c(.1,.9),
        legend.title = element_blank(), 
        legend.background = element_blank()) +
  labs(x="Distance to nearest adaptive locus", y="Proportion of genes")
ggsave("figures/plot_dist_to_goebl_loci_95_quant.pdf", plot = plot_dist_to_goebl_loci_95_quant,
       device = "pdf", height = 65.625, width = 87.5, dpi = 300, units = "mm")

plot_dist_to_goebl_loci_99_quant <- ggplot(data = dist_to_adaptive) +
    geom_line(aes(x=distance, y=proportion, linetype=set, color=set),
              linewidth=.25) +
    scale_color_manual(values = c("red", "black")) +
    scale_linetype_manual(values = c(1,2)) +
    geom_ribbon(data=random_expectation_95CI,
                aes(x=distance, ymin = min, ymax = max),
                color='grey50', alpha=.25, linewidth=.25) +
    scale_x_continuous(limits = c(0,1000001)) +
    scale_y_continuous(limits = c(0,.1)) +
    theme_bw(base_size = 12) +
    theme(legend.position =  c(.1,.9),
          legend.title = element_blank(), 
          legend.background = element_blank()) +
    labs(x="Distance to nearest adaptive locus", y="Proportion of genes")
ggsave("figures/plot_dist_to_goebl_loci_99_quant.pdf", plot = plot_dist_to_goebl_loci_99,
       device = "pdf", height = 65.625, width = 87.5, dpi = 300, units = "mm")

#### randomization test ####
test_control_df <- data.frame(test_ranges) %>% 
  mutate(set="test") %>%
  full_join(data.frame(control_ranges) %>%
              mutate(set="control"))
reps <- 1000
r <- rep(NA,1000)
for (i in 1:length(r) ){
  r[i] <- mean(sample(control_ranges$dist_to_adaptive, 1038))
}
# empirical p-value (one-tailed) is the proportion of randomized means 
# less than or equal to the observed mean
length(r[r<=mean(test_ranges$dist_to_adaptive)]) / length(r)
mean(test_ranges$dist_to_adaptive)
mean(control_ranges$dist_to_adaptive)
