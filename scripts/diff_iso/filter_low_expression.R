#### filter RSEM output to remove genes/isoforms with low expected counts ####

# list of Trinity isoforms matched to Trinity genes
iso_gene_map <- read.table("data2/transcriptome/all_sample_Trinity.cd-hit-est_99.iso_gene_map.txt",
                           col.names = c("gene", "isoform"))

gene_counts_raw <- read.table("data2/RSEM_out_23-2-1/gene_counts_matrix.txt",
                               header = T, row.names = 1)

isoform_counts_raw <- read.table("data2/RSEM_out/isoform_counts_matrix.txt",
                                 header=T, row.names = 1)

# filter at gene level
gene_counts_filtered <- gene_counts_raw %>% 
  mutate(dune_total_counts=rowSums(.[1:12]),
         nondune_total_counts=rowSums(.[13:24]),
         dune_n=rowSums(.[1:12]>=3),
         nondune_n=rowSums(.[13:24]>=3)) %>%
  filter(dune_total_counts >= 24 &
           nondune_total_counts >= 24 &
           dune_n >= 8 & nondune_n >= 8)

low_expressed_genes <- data.frame(gene=rownames(subset(
  gene_counts_raw, !(rownames(gene_counts_raw) %in%
                       rownames(gene_counts_filtered)))))

# filter at isoform level
isoforms_keep <- rowSums(isoform_counts_raw) >= 24

isoform_counts_filtered <- isoform_counts_raw[isoforms_keep,]
dim(isoform_counts_filtered)

low_expressed_isoforms <- data.frame(isoform=rownames(subset(isoform_counts_raw,
                                 !(rownames(isoform_counts_raw) %in%
                                     rownames(isoform_counts_filtered)))))

# combine into single list
low_expr <- left_join(low_expressed_genes, iso_gene_map) %>%
  dplyr::select(isoform) %>%
  bind_rows(low_expressed_isoforms) %>%
  unique()


write.table(low_expr, file = "analysis/diff_iso_out/2023-02-01/low_expr_isoforms.txt",
            quote = F, row.names = F, col.names = F)


