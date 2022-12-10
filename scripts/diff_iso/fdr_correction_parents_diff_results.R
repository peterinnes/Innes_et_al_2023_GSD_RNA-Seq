#### read-in results from parents_diff_v2.py ####
library(dplyr)
diff_iso_results <- read.table("analysis/diff_iso_out/parents_diff_iso_results.2022-10-10.txt",
                               sep = ' ', col.names = c("Trinity_gene", "pval"))
diff_iso_results$fdr <- p.adjust(diff_iso_results$pval, method="fdr",
                                 n = nrow(diff_iso_results)) 

ds_genes <- subset(diff_iso_results, fdr<.05) %>%
  arrange(fdr)

write.table(ds_genes,
            "analysis/diff_iso_out/signigicant_parents_diff_iso_results.2022-10-10.txt",
            sep='\t', quote = F, row.names = F)

subset(isoform_counts_raw, rownames(isoform_counts_raw) %in% c("TRINITY_DN90709_c13_g1_i2","TRINITY_DN90709_c13_g1_i1"))

subset(isoform_counts_raw, rownames(isoform_counts_raw) %in% c("TRINITY_DN1301_c0_g1_i10","TRINITY_DN1301_c0_g1_i6"))
