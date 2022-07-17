#### analysis of DE/DS genes in inversion regions ####
library(ggplo2)
library(ggstatsplot)

# read in table of inversion positions/info (Table 1 from Huang et al 2020)
invr_regions <- read.table("inversion_regions.txt", col.names = c("MDS", "chr", "start", "end", "num_of_outlier_wind", "PC1_varianace", "PC2_variance", "proportion_of_between_cluster_sum_of_squares", "region_code")) %>%
  dplyr::select(region_code, chr, start, end)


# Read in the annotations and Filter for only expressed genes)
gff <- data.frame(rtracklayer::import("data/ref_genome_Ha412HO/HAN412_Eugene_curated_v1_1.gff3")) %>%
  rename(seqnames="chrom", ID="Ha412_gene")

expressed_genes_gff <- subset(gff, type=="gene")[1:11] %>% inner_join(expressed_genes)

write.table(expressed_genes_gff, file="data/ref_genome_Ha412HO/HAN412_Eugene_curated_v1_1.expressed_genes_gff.tmp", sep = "\t", row.names = F, col.names = F, quote = F)



# maybe first use bedtools map -o count_distinct to count the number of genes from gff file in each inversion?
# https://www.biostars.org/p/65195/

# need to subtract 1 from start position when converting from gff to bed bc gff is 1-based coordinate system while bed is 0-based
DE_gene_regions.bed <- left_join(de_results_p05_Shrink, expressed_genes_gff) %>%
  dplyr::select(chrom, start, end) %>%
  mutate(start=start-1)

write.table(DE_gene_regions.bed, file = "analysis/DESeq2/DE_gene_regions.bed", sep = "\t", quote = F, row.names = F, col.names = F)

#### use bedtools and basic maths to get counts of DE genes within inversions, non-DE genes within inversions, DE genes outside inversions, non-DE genes outside inversions 

# count DE genes within inversions = 760 
# bedtools intersect -loj -a analysis/inversions/inversion_regions.bed -b analysis/DESeq2/DE_gene_regions.bed | wc -l

# count number of expressed genes of any kind in the inversions = 3034; non-DE genes within inversions = 3034 - 760 = 2274
# bedtools map -a analysis/inversions/inversion_regions.bed -b data/ref_genome_Ha412HO/HAN412_Eugene_curated_v1_1.expressed_genes_gff.tmp -c 10 -o count_distinct | cut -f4 | paste -sd+ | bc

# count number of DE genes outside inversions = 4320
# bedtools intersect -v -a analysis/DESeq2/DE_gene_regions.bed -b analysis/inversions/inversion_regions.bed | wc -l

# count number of expressed genes of any kind outside inversions = 29086; non-DE genes outside inversions = 29086 - 4320 = 24766
# bedtools intersect -v -a data/ref_genome_Ha412HO/HAN412_Eugene_curated_v1_1.expressed_genes_gff.tmp -b analysis/inversions/inversion_regions.bed | grep -v 'Chr00' | cut -f1 | sort | uniq -c | tr -s ' ' | cut -d ' ' -f 2 | paste -sd+ | bc

DE_inv_enrichment_table <- data.frame("DE" = c(760, 4320), "non-DE" = c(3034,24766),
                               row.names = c("inversion", "non-inversion"),
                               stringsAsFactors = F)

DE_inv_fe_test <- fisher.test(DE_inv_enrichment_table)

x <- c()
for (row in rownames(DE_inv_enrichment_table)) {
  for (col in colnames(DE_inv_enrichment_table)) {
    x <- rbind(x, matrix(rep(c(row, col), DE_inv_enrichment_table[row, col]), ncol = 2, byrow = TRUE))
  }
}
DE_inv_enrichment_df <- as.data.frame(x)
colnames(DE_inv_enrichment_df) <- c("Inversion", "DE")

ggbarstats(
  DE_inv_enrichment_df, DE, Inversion,
  results.subtitle = FALSE,
  subtitle = paste0(
    "Fisher's exact test", ", p-value = ",
    ifelse(DE_inv_fe_test$p.value < 0.001, "< 0.001", round(DE_inv_fe_test$p.value, 3))
  )
)

#### Enrichment of DS genes? ####

# DEXSeq (DEU)
deu_gene_regions.bed <- left_join(deu_genes_noLFCthreshold, expressed_genes_gff) %>%
  dplyr::select(chrom, start, end) %>%
  mutate(start=start-1)

write.table(deu_gene_regions.bed, file = "analysis/DEXSeq/deu_gene_regions.bed", sep = "\t", quote = F, row.names = F, col.names = F)

# bash commands

# Count DEU  genes within inversions = 273
# bedtools intersect -loj -a ~/gsd_RNA-seq/analysis/inversions/inversion_regions.bed -b analysis/DEXSeq/deu_gene_regions.bed | wc -l

# count number of genes of any kind in the inversions = 3034 ; non-DEU genes within inversions = 3034 - 273 = 2761
# bedtools map -a analysis/inversions/inversion_regions.bed -b data/ref_genome_Ha412HO/HAN412_Eugene_curated_v1_1.expressed_genes_gff.tmp -c 10 -o count_distinct | cut -f4 | paste -sd+ | bc

# count number of DEU genes outside inversions = 683
# bedtools intersect -v -a analysis/DEXSeq/deu_gene_regions.bed -b analysis/inversions/inversion_regions.bed | wc -l

#  count number of genes of any kind outside inversions = 42056; non-DEU genes outside inversions = 29086 - 683 = 28403
# bedtools intersect -v -a data/ref_genome_Ha412HO/HAN412_Eugene_curated_v1_1.expressed_genes_gff.tmp -b analysis/inversions/inversion_regions.bed | grep -v 'Chr00' | cut -f1 | sort | uniq -c | tr -s ' ' | cut -d ' ' -f 2 | paste -sd+ | bc

deu_inv_enrichment_table <- data.frame("DEU" = c(273, 683), "non-DEU" = c(2761,28403),
                               row.names = c("inversion", "non-inversion"),
                               stringsAsFactors = F)

deu_inv_fe_test <- fisher.test(deu_inv_enrichment_table)

# rMATS
rmats_ds_gene_regions.bed <- left_join(rmats_ds_genes, expressed_genes_gff) %>%
  dplyr::select(chrom, start, end) %>%
  mutate(start=start-1)
write.table(rmats_ds_gene_regions.bed, file = "analysis/rMATS/results_2022-07-14/rmats_ds_gene_regions.bed", sep = "\t", quote = F, row.names = F, col.names = F)

# number of DS genes within inversions
#bedtools intersect -loj -a ~/gsd_RNA-seq/analysis/inversions/inversion_regions.bed -b analysis/rMATS/results_2022-07-14/rmats_ds_gene_regions.bed | wc -l
#159

# count number of expressed genes of any kind in the inversions = 3034 ; non-DS genes within inversions = 3034 - 159 = 2875
# bedtools map -a analysis/inversions/inversion_regions.bed -b data/ref_genome_Ha412HO/HAN412_Eugene_curated_v1_1.expressed_genes_gff.tmp -c 10 -o count_distinct | cut -f4 | paste -sd+ | bc

# number ds genes outside of inversions = tot ds genes - ds genes inside inversions = 1038 - 159 = 879

#  count number of expressed genes of any kind outside inversions = 42056; non-DS genes outside inversions = 29086 - 879 = 28207
# bedtools intersect -v -a data/ref_genome_Ha412HO/HAN412_Eugene_curated_v1_1.expressed_genes_gff.tmp -b analysis/inversions/inversion_regions.bed | grep -v 'Chr00' | cut -f1 | sort | uniq -c | tr -s ' ' | cut -d ' ' -f 2 | paste -sd+ | bc

ds_inv_enrichment_table <- data.frame("DS" = c(159, 879), "non-DS" = c(2875,28207),
                                       row.names = c("inversion", "non-inversion"),
                                       stringsAsFactors = F)

ds_inv_fe_test <- fisher.test(ds_inv_enrichment_table)
