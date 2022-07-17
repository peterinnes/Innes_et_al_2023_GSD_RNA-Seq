# eQTL/sQTL analysis with Matrix eQTL

library(MatrixEQTL)

#### get input data ####
snp_mat <- read.table("data/hc_out/dune_non-dune.nHet_filtered.hc.012", header = F, row.names = 1) # read-in ¨012¨ format snp matrix. takes a while bc R things
indvs <- read.table("data/hc_out/dune_non-dune.nHet_filtered.hc.012.indv", col.names = "sample") #sample names
snps_pos <- read.table("data/hc_out/dune_non-dune.nHet_filtered.hc.012.pos", col.names = c("chr", "pos")) %>%
  rownames_to_column("SNP") %>% #snp ids and positions
  mutate(SNP=as.integer(SNP))
# get positions of the genes. using 'transcripts' object from analyze_differential_expression_DESeq2.R.
genes_pos <- dplyr::select(transcripts, ID, chrom, start, end) #genes positions

# get the matrix sorted out with sample names and snp ids, and transpose so samples are columns
row.names(snp_mat) <- indvs$sample
snp_mat <- data.frame(t(snp_mat))
rownames(snp_mat) <- snps_pos$SNP
snp_mat <- as.matrix(snp_mat)

# convert snp_mat into matrixEQTL 'SlicedData' format
snps <- SlicedData$new()
snps$initialize(snp_mat)

# get expression data from analyze_differential_expression_DESeq2.R and convert to matrixEQTL format.
expr_mat <- as.matrix(normalized_counts)
genes <- SlicedData$new()
genes$initialize(expr_mat)

# get covariate data (PC1 from SNP PCA) from pca.R script and convert to matrixEQTL format
pc_mat <- dplyr::select(snp_pca_df, EV1)
rownames(pc_mat) <- snp_pca_df$sample.id
pc_mat <- t(as.matrix(pc_mat))
cvrt <- SlicedData$new()
cvrt$initialize(pc_mat)

#### run eQTL ####

me <- Matrix_eQTL_main(
  snps = snps,
  gene = genes,
  cvrt = cvrt,
  output_file_name = "analysis/matrixEQTL/eQTL_trans.out",
  pvOutputThreshold = 0, # set trans threshold to zero to only consider cis eQTL
  useModel = modelLINEAR,
  errorCovariance = numeric(),
  verbose = TRUE,
  output_file_name.cis = "analysis/matrixEQTL/eQTL_cis.out",
  pvOutputThreshold.cis = 2e-2,
  snpspos = snps_pos,
  genepos = genes_pos,
  cisDist = 1e6,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE)

# "Matching data files and location files
# 32308 of 32311 genes matched
# 295383 of 295383 SNPs matched
# Results
cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
cat('Detected local eQTLs:', '\n');
show(me$cis$eqtls)
cat('Detected distant eQTLs:', '\n');
show(me$trans$eqtls)

# check the qqplot
plot(me)

# read-in cis-eQTL results
cis_eQTL <- read.table("analysis/matrixEQTL/eQTL_cis.out", header = T) %>%
  left_join(snps_pos)
cis_eQTL_sig <- subset(cis_eQTL, FDR < .05)

# 8,267 SNPs associated with expression of 2578 genes
dim(cis_eQTL_sig)
length(unique(cis_eQTL_sig$gene))

trans_eQTL <- read.table("analysis/matrixEQTL/eQTL_trans.out", header = T) %>%
  left_join(snps_pos) %>%
  subset(FDR < .05)

#### splicing QTL ####

# get splicing "phenotype" data, from rMATS
