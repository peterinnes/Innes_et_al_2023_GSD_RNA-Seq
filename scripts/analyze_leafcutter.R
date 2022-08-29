BiocManager::install("DirichletMultinomial")
devtools::install_github("davidaknowles/leafcutter/leafcutter")
install.packages("rstan")
install.packages("rstantools")

library(rstan)
library(leafcutter)


# read in intron excision ratio data. These data have been normalized,
# via the leafcutter script prepare_phenotype_table.py
ier_df <- read.table("data2/leafcutter_out/dune_vs_non-dune_perind.counts.gz.qqnorm_chr1", header = T)

# results are split into separate files per chromosome, so we need to combine them
for(i in 2:17){
  ier_df <- dplyr::bind_rows(ier_df, read.table(paste0("data2/leafcutter_out/dune_vs_non-dune_perind.counts.gz.qqnorm_chr",i), header = T))
}

rownames(ier_df) <- ier_df$ID
