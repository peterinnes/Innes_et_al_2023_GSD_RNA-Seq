#### script for filtering and normalizing/transforming raw counts ####
#Written by Peter Innes 

library(DESeq2)
library(dplyr)
library(matrixStats)

sample_table_deseq2 <- read.table("data/expData/sample_table.txt", header=T)
sample_table_deseq2$habitat <- as.factor(sample_table_deseq2$habitat)
counts <- read.table("data/expData/raw_counts.txt",
                     row.names = 1)

names(counts) <- paste(sample_table_deseq2$habitat,"_",
                       sample_table_deseq2$sample_no, sep = "")


# pre-Filter out lowly expressed genes,
# where read count is less than 24 across all samples combined
# This is were you could filter by median or mean read count instead, 
# I think using rowMeans() or something
#keep <- rowSums(counts(dds)) >= 24
keep <- which(rowMeans(counts) >= 1) #same as line above but doesn't require dds obj
#Here I filter out genes with expression counts with a median less than 5 
filtered_counts <- counts[keep,]

# do the regularized-log transformation, used in expression PCA
# and also used in  WGCNA potentially. Needs to be applied
# WITHOUT normalization bc it requires expression values to be integers
rlog_filtered_counts <- rlog(as.matrix(filtered_counts)) %>%
  data.frame() %>%
  tibble::rownames_to_column("Ha412_gene")

# write to file
write.table(rlog_filtered_counts, file="../data/Rdata/rlog-transformed_expression_counts.tsv",
            quote = F, row.names = F, sep = "\t")
