#### read-in results from parents_diff_v2.py ####
library(dplyr)
diff_iso_results <- read.table("analysis/diff_iso_out/2023-02-01/parents_diff_iso_results.2023-02-04.txt",
                               sep = ' ', col.names = c("Trinity_gene", "pval"))
diff_iso_results$fdr <- p.adjust(diff_iso_results$pval, method="fdr",
                                 n = nrow(diff_iso_results)) 

ds_genes <- subset(diff_iso_results, fdr<.05) %>%
  arrange(fdr)

write.table(ds_genes,
            "analysis/diff_iso_out/2023-02-01/significant_parents_diff_iso_results.2023-02-05.txt",
            sep='\t', quote = F, row.names = F)
