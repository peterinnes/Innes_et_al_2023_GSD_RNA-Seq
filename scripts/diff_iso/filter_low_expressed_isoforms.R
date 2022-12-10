isoform_counts_raw <- read.table("data2/RSEM_out/isoform_counts_matrix.txt",
                                 header=T, row.names = 1)

isoforms_keep <- rowSums(isoform_counts_raw) >= 24

isoform_counts_filtered <- isoform_counts_raw[isoforms_keep,]
dim(isoform_counts_filtered)

low_expressed_isoforms <- subset(isoform_counts_raw,
                                 !(rownames(isoform_counts_raw) %in%
                                     rownames(isoform_counts_filtered)))

write.table(rownames(low_expressed_isoforms),
            file = "analysis/diff_iso_out/low_expr_isoforms.txt", quote = F,
            row.names = F)


