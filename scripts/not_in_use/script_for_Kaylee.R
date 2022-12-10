#### script for filtering and normalizing/transforming raw counts ####

library(DESeq2)
library(dplyr)

sample_table_deseq2 <- read.table("analysis/DESeq2/sample_table.txt", header=T)
sample_table_deseq2$habitat <- as.factor(sample_table_deseq2$habitat)
counts <- read.table("analysis/DESeq2/htseq-count_out/htseq-count_results.2022-6-26.txt",
                     row.names = 1)

names(counts) <- paste(sample_table_deseq2$habitat,"_",
                       sample_table_deseq2$sample_no, sep = "")


# pre-Filter out lowly expressed genes,
# where read count is less than 24 across all samples combined
# This is were you could filter by median or mean read count instead, 
# I think using rowMeans() or something
keep <- rowSums(counts) >= 24
filtered_counts <- counts[keep,]

# do the regularized-log transformation, used in expression PCA
# and alsoused in  WGCNA potentially. Needs to be applied
# WITHOUT normalization bc it requires expression values to be integers
rlog_filtered_counts <- rlog(as.matrix(filtered_counts)) %>%
  data.frame() %>%
  tibble::rownames_to_column("Ha412_gene")

# write to file
write.table(rlog_filtered_counts, file="analysis/DESeq2/rlog-transformed_expression_counts.tsv",
            quote = F, row.names = F, sep = "\t")
