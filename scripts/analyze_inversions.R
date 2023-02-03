######## analysis of DE/DS genes in inversion regions
library(ggplot2)
library(dplyr)
library(ggstatsplot)
library(patchwork)

#### Read in data ####
# read in table of inversion positions/info (Table 1 from Huang et al 2020)
inv_regions <- read.table("inversion_regions.txt",
                          col.names = c("MDS", "chr", "start", "end",
                                        "num_of_outlier_wind", "PC1_varianace",
                                        "PC2_variance",
                                        "proportion_of_between_cluster_sum_of_squares",
                                        "region_code")) %>%
  dplyr::select(region_code, chr, start, end)


# Read in the annotations and filter for only expressed genes)
gff <- data.frame(rtracklayer::import("data/ref_genome_Ha412HO/HAN412_Eugene_curated_v1_1.gff3")) %>%
  rename(chrom="seqnames", Ha412_gene="ID")

expressed_genes <- read.table("analysis/DESeq2/expressed_genes.txt",
                              col.names = "Ha412_gene")

expressed_genes_gff <- subset(gff, type=="gene")[1:11] %>% inner_join(expressed_genes)

write.table(expressed_genes_gff,
            file="data/ref_genome_Ha412HO/HAN412_Eugene_curated_v1_1.expressed_genes_gff.tmp",
            sep = "\t", row.names = F, col.names = F, quote = F)


#### Enrichment of DE genes within inversions ####
DE_genes <- read.table("analysis/GO_analysis/study_DE_genes_noLFCthreshold.txt",
                       col.names = "Ha412_gene")

# need to subtract 1 from start position when converting from gff to bed bc gff is 1-based coordinate system while bed is 0-based
DE_gene_regions.bed <- left_join(DE_genes, expressed_genes_gff) %>%
  dplyr::select(chrom, start, end) %>%
  mutate(start=start-1)

write.table(DE_gene_regions.bed, file = "analysis/DESeq2/DE_gene_regions.bed",
            sep = "\t", quote = F, row.names = F, col.names = F)

# use bedtools and basic maths to get counts of DE genes within inversions, non-DE genes within inversions, DE genes outside inversions, non-DE genes outside inversions 

# count DE genes within inversions = 593
# bedtools intersect -loj -a analysis/inversions/four_inversion_regions.bed -b analysis/DESeq2/DE_gene_regions.bed | wc -l

# count number of expressed genes of any kind in the inversions = 1958; non-DE genes within inversions = 1958 - 593 = 1365
# bedtools map -a ~/gsd_RNA-seq/analysis/inversions/four_inversion_regions.bed -b ~/gsd_RNA-seq/data/ref_genome_Ha412HO/HAN412_Eugene_curated_v1_1.expressed_genes_gff.tmp -c 10 -o count_distinct | cut -f4 | paste -sd+ | bc

# count number of DE genes outside inversions = 4510
# bedtools intersect -v -a ~/gsd_RNA-seq/analysis/DESeq2/DE_gene_regions.bed -b ~/gsd_RNA-seq/analysis/inversions/four_inversion_regions.bed | wc -l

# count number of expressed genes of any kind outside inversions = 30350; non-DE genes outside inversions = 30350 - 4510 = 25840
# bedtools intersect -v -a ~/gsd_RNA-seq/data/ref_genome_Ha412HO/HAN412_Eugene_curated_v1_1.expressed_genes_gff.tmp -b ~/gsd_RNA-seq/analysis/inversions/four_inversion_regions.bed | cut -f1 | sort | uniq -c | tr -s ' ' | cut -d ' ' -f 2 | paste -sd+ | bc

DE_inv_enrichment_table <- data.frame("DE" = c(593, 4510), "non-DE" = c(1365,25840),
                               row.names = c("Inversion", "Non-inversion"),
                               stringsAsFactors = F)

DE_inv_fe_test <- fisher.test(DE_inv_enrichment_table)

x <- c()
for (row in rownames(DE_inv_enrichment_table)) {
  for (col in colnames(DE_inv_enrichment_table)) {
    x <- rbind(x, matrix(rep(c(row, col), DE_inv_enrichment_table[row, col]), ncol = 2, byrow = TRUE))
  }
}
DE_inv_enrichment_df <- as.data.frame(x) %>%
  rename(Inversion_status="V1", Expression_status="V2") %>%
  mutate(Expression_status=gsub("\\.","-",Expression_status))

DE_inv_enrichment_df$Inversion_status <- as.factor(DE_inv_enrichment_df$Inversion_status)
DE_inv_enrichment_df$Expression_status <- as.factor(DE_inv_enrichment_df$Expression_status)

plot_inv_enrich_DE <- ggplot(DE_inv_enrichment_df,
       aes(fill=factor(Expression_status,
                       levels=c("non-DE", "DE")),
           x=Inversion_status)) + 
  geom_bar(position="fill", alpha=0.5) +
  scale_fill_discrete(breaks=c("DE", "non-DE"), type = c("grey50","#317EC2")) +
  theme_classic(base_size = 12) +
  theme(legend.title = element_blank(),
        legend.position = "top") +
  labs(x="", y="Proportion")#, title = "Fisher's exact test, p < .001")


#### Enrichment of DS genes within inversion ####
# rMATS
rmats_ds_gene_regions.bed <- left_join(rmats_ds_genes, expressed_genes_gff) %>%
  dplyr::select(chrom, start, end) %>%
  mutate(start=start-1)

write.table(rmats_ds_gene_regions.bed, file = "analysis/rMATS/results_2022-07-14/rmats_ds_gene_regions.bed", sep = "\t", quote = F, row.names = F, col.names = F)

# number of DS genes within inversions = 129
# bedtools intersect -loj -a ~/gsd_RNA-seq/analysis/inversions/four_inversion_regions.bed -b ~/gsd_RNA-seq/analysis/rMATS/results_2022-07-14/rmats_ds_gene_regions.bed | wc -l

# count number of expressed genes of any kind in the inversions = 1958 ; non-DS genes within inversions = 1958 - 129 = 1829
# bedtools map -a ~/gsd_RNA-seq/analysis/inversions/four_inversion_regions.bed -b ~/gsd_RNA-seq/data/ref_genome_Ha412HO/HAN412_Eugene_curated_v1_1.expressed_genes_gff.tmp -c 10 -o count_distinct | cut -f4 | paste -sd+ | bc

# count number of DS genes outside of inversions = tot ds genes - ds genes inside inversions = 1038 - 129 = 909

#  count number of expressed genes of any kind outside inversions = 30350; count non-DS genes outside inversions = 30350 - 909 = 29441
# bedtools intersect -v -a ~/gsd_RNA-seq/data/ref_genome_Ha412HO/HAN412_Eugene_curated_v1_1.expressed_genes_gff.tmp -b ~/gsd_RNA-seq/analysis/inversions/inversion_regions.bed | cut -f1 | sort | uniq -c | tr -s ' ' | cut -d ' ' -f 2 | paste -sd+ | bc

DS_inv_enrichment_table <- data.frame("DS" = c(129, 909), "non-DS" = c(1829,29441),
                                      row.names = c("Inversion", "Non-inversion"),
                                      stringsAsFactors = F)

DS_inv_fe_test <- fisher.test(DS_inv_enrichment_table)

x <- c()
for (row in rownames(DS_inv_enrichment_table)) {
  for (col in colnames(DS_inv_enrichment_table)) {
    x <- rbind(x, matrix(rep(c(row, col),
                             DS_inv_enrichment_table[row, col]), ncol = 2, byrow = TRUE))
  }
}

DS_inv_enrichment_df <- as.data.frame(x) %>%
  rename(Inversion_status="V1", Splicing_status="V2") %>%
  mutate(Splicing_status=gsub("\\.","-",Splicing_status))

DS_inv_enrichment_df$Inversion_status <- as.factor(DS_inv_enrichment_df$Inversion_status)
DS_inv_enrichment_df$Splicing_status <- as.factor(DS_inv_enrichment_df$Splicing_status)

plot_inv_enrich_DS <- ggplot(DS_inv_enrichment_df,
       aes(fill=factor(Splicing_status,
                       levels=c("non-DS", "DS")),
           x=Inversion_status)) + 
  geom_bar(position="fill", alpha=0.5) +
  scale_fill_discrete(breaks=c("DS", "non-DS"), type = c("grey50","red")) +
  theme_classic(base_size = 12) +
  theme(legend.title = element_blank(),
        legend.position = "top") +
  labs(x="", y="")#, y="Proportion")#, title = "Fisher's exact test, p < .001")

plot_inv_enrich <- plot_inv_enrich_DE | plot_inv_enrich_DS

ggsave("figures/plot_inversion_enrich_raw.pdf", plot=plot_inv_enrich,
       device = "pdf", height = 87.5 , width = 131.25, dpi = 300, units = "mm")


#### look at DS and DE genes within seed size qtl ####
system("bedtools intersect -loj -a ~/gsd_RNA-seq/seed_size_qtl_pos.txt -b ~/gsd_RNA-seq/analysis/rMATS/results_2022-07-14/rmats_ds_gene_regions.bed > ~/gsd_RNA-seq/analysis/rMATS/results_2022-07-14/ds_genes_within_seed_size_qtl.txt")

ds_seed_size_qtl_genes <- read.table("analysis/rMATS/results_2022-07-14/ds_genes_within_seed_size_qtl.txt") %>%
  rename(qtl_chrom=V1, qtl_start=V2, qtl_end=V3, chrom=V4, start=V5, end=V6) %>%
  mutate(start=start+1) %>% #add one to convert from bed to gff
  left_join(expressed_genes_gff) %>%
  dplyr::select(Ha412_gene)

# DS_embryo_dev_genes from analyze_GO_enrichment.R
filter(all_AS_events_deltaPSI, GeneID %in% DS_embryo_dev_genes) %>%
  filter(GeneID %in% ds_seed_size_qtl_genes$Ha412_gene) %>%
  filter(FDR <.05 )

slice_max(filter(all_AS_events_deltaPSI, GeneID %in% ds_seed_size_genes$Ha412_gene) %>%
  filter(FDR <.05), IncLevelDifference, n=5)

slice_min(filter(all_AS_events_deltaPSI, GeneID %in% ds_seed_size_genes$Ha412_gene) %>%
            filter(FDR <.05), IncLevelDifference, n=5)

system("bedtools intersect -loj -a ~/gsd_RNA-seq/seed_size_qtl_pos.txt -b ~/gsd_RNA-seq/analysis/DESeq2/DE_gene_regions.bed > ~/gsd_RNA-seq/analysis/DESeq2/DE_genes_within_seed_size_qtl.txt")
de_seed_size_qtl_genes <- read.table("analysis/DESeq2/DE_genes_within_seed_size_qtl.txt") %>%
  rename(qtl_chrom=V1, qtl_start=V2, qtl_end=V3, chrom=V4, start=V5, end=V6) %>%
  mutate(start=start+1) %>% #add one to convert from bed to gff
  left_join(expressed_genes_gff) %>%
  dplyr::select(Ha412_gene)

# DE_embryo_seed_dev_genes from analyze_GO_enrichment.R
filter(de_results_Shrink_df, Ha412HOv2_gene %in% DE_embryo_seed_dev_genes) %>%
  filter(Ha412HOv2_gene %in% de_seed_size_qtl_genes$Ha412_gene)
