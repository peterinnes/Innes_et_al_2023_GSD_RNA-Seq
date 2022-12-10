# eQTL/sQTL analysis with Matrix eQTL

library(MatrixEQTL)

#### get input data ####
snp_mat <- read.table("data/hc_out/dune_non-dune.nHet_filtered.hc.012",
                      header = F, row.names = 1) #load ¨012¨ format snp matrix.
indvs <- read.table("data/hc_out/dune_non-dune.nHet_filtered.hc.012.indv",
                    col.names = "sample") #sample names
snps_pos <- read.table("data/hc_out/dune_non-dune.nHet_filtered.hc.012.pos",
                       col.names = c("chr", "pos")) %>%
  rownames_to_column("SNP") %>% #snp ids and positions
  mutate(SNP=as.integer(SNP))
# get positions of the genes. using 'transcripts' object from analyze_differential_expression_DESeq2.R.

genes_pos <- dplyr::select(transcripts, ID, chrom, start, end) #genes positions

# get the matrix sorted out with sample names and snp ids,
# and transpose so samples are columns
row.names(snp_mat) <- indvs$sample
snp_mat <- data.frame(t(snp_mat))
rownames(snp_mat) <- snps_pos$SNP
snp_mat <- as.matrix(snp_mat)

# convert snp_mat into matrixEQTL 'SlicedData' format
snps <- SlicedData$new()
snps$initialize(snp_mat)

# get expression data from analyze_differential_expression_DESeq2.R
# and convert to matrixEQTL format.
expr_mat <- as.matrix(normalized_counts)
genes <- SlicedData$new()
genes$initialize(expr_mat)

# get covariate data (PC1 from SNP PCA) from pca.R script
# and convert to matrixEQTL format
pc_mat <- dplyr::select(snp_pca_df, EV1)
rownames(pc_mat) <- snp_pca_df$sample.id
pc_mat <- t(as.matrix(pc_mat))
cvrt <- SlicedData$new()
cvrt$initialize(pc_mat)

#### expression QTL ####
me <- Matrix_eQTL_main(
  snps = snps,
  gene = genes,
  cvrt = cvrt,
  output_file_name = "analysis/matrixEQTL/eQTL_trans.out",
  pvOutputThreshold = 0, #set trans threshold to zero to only consider cis eQTL
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
  left_join(snps_pos) %>%
  rename(chr="chrom")
cis_eQTL_sig <- filter(cis_eQTL, FDR < .05)

# 8,267 SNPs associated with expression of 2578 genes
dim(cis_eQTL_sig)
length(unique(cis_eQTL_sig$gene))

-log10(max(cis_eQTL_sig$p.value)) #get threshold for significance in pvalue scale

#trans_eQTL <- read.table("analysis/matrixEQTL/eQTL_trans.out", header = T) %>%
#  left_join(snps_pos) %>%
#  subset(FDR < .05)

#### eQTL manhattan plot ####
chrom_lengths <- read.table("data/ref_genome_Ha412HO/chrom_sizes_Ha412HOv2.0.txt",
                            col.names = c("chrom", "length")) %>%
  mutate(length=as.numeric(length))

cum_lengths <- chrom_lengths %>%
  mutate(cumstart=cumsum(length)-length,
         cumend=cumsum(length))

# get inversion regions
inversion_regions <- read.table("analysis/inversions/inversion_regions.txt") %>%
  dplyr::select(V2,V3,V4,V9) %>%
  rename("chrom"=V2, "start"=V3, "end"=V4, "name"=V9) %>%
  filter(name %in% c("pet05.01", "pet09.01", "pet11.01", "pet17.01")) %>%
  left_join(cum_lengths[c(1,3)]) %>%
  mutate(inv_cumstart=start+cumstart, inv_cumend=end+cumstart)

# randomly sample just 50% of non-sig SNPs, for improved visualization
cis_eQTL_reduced <- filter(cis_eQTL, FDR >= .05) %>%
  sample_frac(.5) %>%
  bind_rows(filter(cis_eQTL, FDR <.05)) %>%
  inner_join(cum_lengths) %>%
  mutate(cum_pos = pos + cumstart)

axisdf_eQTL <- cis_eQTL_reduced %>%
  group_by(chrom) %>%
  summarize(center = (max(cum_pos) + min(cum_pos)) / 2) %>%
  mutate(chrom=c("1","2","3","4","5",
                 "6","7","8","9","10",
                 "11","12","13","14","15",
                 "16","17"))

eQTL_manhattan_plot <- ggplot(cis_eQTL_reduced,
                              aes(x=cum_pos, y=-log10(p.value),
                                  color=as.factor(chrom))) +
  
  ## mark inversion regions
  #geom_rect(data=inversion_regions,
  #          aes(xmin=inv_cumstart, xmax=inv_cumend,
  #              ymin=-1,
  #              ymax=0),
  #          fill="red", inherit.aes = F) +
  
  # add each SNP–gene correlation.
  geom_point(shape=1) +
  scale_color_manual(values = rep(c("grey50", "black"),
                                  unique(length(axisdf_eQTL$chrom)))) +
  #scale_alpha_manual(name="FDR < .05",values = c(.1,1), labels=c("no", "yes")) +
  #scale_size_manual(name="FDR < .05", values = c(1,4), labels=c("no", "yes")) +
  #scale_color_gradientn(colors=c("grey","turquoise","blue","navyblue","black")) +
  
  # custom axes, adjust limits of y axis.
  scale_x_continuous(label = axisdf_eQTL$chrom,
                     breaks = axisdf_eQTL$center, expand=c(.01,.01)) +
  scale_y_continuous(expand=c(.01,.01)) +
  ylim(1.7, max(-log10(cis_eQTL_reduced$p.value))) +
  coord_cartesian(clip = 'off') +
  
  geom_hline(lty=2, color="red", yintercept=4.35) +
  
  # custom theme
  theme_bw() +
  theme(legend.position = "none",
        text = element_text(size=24),
        axis.line.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x=element_text(size=18),
        plot.title = element_text(size=18),
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  labs(x="", y="-log10 p-value", title="cis-eQTL")

eQTL_manhattan_plot

#### splicing QTL ####

# get splice data 

rmats_splice_pheno_df <- rmats_splice_pca_df
colnames(rmats_splice_pheno_df) <- colnames(rmats_splice_pheno_df) %>%
    #gsub("_R1_001.trimmed.fq.gz","",.) %>%
    gsub("Kane\\.602\\.","Kane-602-",.)
# reorder columns to be consistent with order in vcf file (same ordering as for eQTL)
rmats_splice_pheno_df <- rmats_splice_pheno_df %>% 
    dplyr::select(as.factor(sample_table_deseq2$fastq))
  
rmats_splice_mat <- as.matrix(rmats_splice_pheno_df)
splice_events <- SlicedData$new()
splice_events$initialize(rmats_splice_mat)



splice_events_pos <- dplyr::select(all_AS_events_deltaPSI, ID, chrom,
                                   exon_start, exon_end) %>%
  subset(ID %in% rownames(rmats_splice_mat))

leafcutter_splice_mat <- as.matrix(ier_df[5:28])

ier_pos <- data.frame(rownames=rownames(leafcutter_splice_mat)) %>%
  separate(col = rownames, sep = ":", into = c("chrom", "start",
                                               "end", "ID")) %>%
  mutate(start=as.numeric(start), end=as.numeric(end)) %>%
  relocate(ID, .before = chrom)

rownames(leafcutter_splice_mat) <- ier_pos$ID
iers <- SlicedData$new()
iers$initialize(leafcutter_splice_mat)

snps_pos_lc <- snps_pos %>%
  mutate(chr=gsub("Ha412HOC","c",chr)) %>%
  mutate(chr=gsub("r0","r",chr))

ms <- Matrix_eQTL_main(
  snps = snps,
  gene = splice_events,
  cvrt = cvrt,
  output_file_name = "analysis/matrixEQTL/sQTL_trans.out",
  pvOutputThreshold = 0, # set trans threshold to zero to only consider cis eQTL
  useModel = modelLINEAR,
  errorCovariance = numeric(),
  verbose = TRUE,
  output_file_name.cis = "analysis/matrixEQTL/sQTL_cis.out",
  pvOutputThreshold.cis = 2e-2,
  snpspos = snps_pos,
  genepos = splice_events_pos,
  cisDist = 1e6,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE)

plot(ms)

# run with leafcutter splicing data
ms_leafcutter <- Matrix_eQTL_main(
  snps = snps,
  gene = iers,
  cvrt = cvrt,
  output_file_name = "analysis/matrixEQTL/lc_sQTL_trans.out",
  pvOutputThreshold = 0, #set trans threshold to zero to only consider cis eQTL
  useModel = modelLINEAR,
  errorCovariance = numeric(),
  verbose = TRUE,
  output_file_name.cis = "analysis/matrixEQTL/lc_sQTL_cis.out",
  pvOutputThreshold.cis = 2e-2,
  snpspos = snps_pos_lc,
  genepos = ier_pos,
  cisDist = 1e6,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE)

# read-in cis-eQTL results
cis_sQTL <- read.table("analysis/matrixEQTL/sQTL_cis.out", header = T) %>%
  rename(gene="ID") %>%
  left_join(splice_events_pos)

cis_sQTL_sig <- filter(cis_sQTL, FDR < .05)

# 11585 SNPs associated with 2747 unique splice events
dim(cis_sQTL_sig)
length(unique(cis_sQTL_sig$ID))

# read-in leafcutter cis-eQTL results
lc_cis_sQTL <- read.table("analysis/matrixEQTL/lc_sQTL_cis.out", header = T) %>%
  rename("ID"=gene) %>%
  left_join(ier_pos)

lc_cis_sQTL_sig <- filter(lc_cis_sQTL, FDR < .05)
dim(lc_cis_sQTL_sig)
length(unique(lc_cis_sQTL_sig$ID))

#### sQTL manhattan plot ####
cis_sQTL <- cis_sQTL %>%
  inner_join(cum_lengths) %>%
  mutate(cum_pos = exon_start + cumstart)

lc_cis_sQTL <- lc_cis_sQTL %>%
  inner_join(cum_lengths %>% 
               mutate(chrom=gsub("Ha412HOC","c",chrom)) %>%
               mutate(chrom=gsub("r0","r",chrom))) %>%
  mutate(chrom=factor(chrom, levels = str_sort(unique(chrom), numeric = T)),
         cum_pos = start + cumstart) %>%
  arrange(chrom)


axisdf_sQTL <- cis_sQTL %>%
  group_by(chrom) %>%
  summarize(center = (max(cum_pos) + min(cum_pos)) / 2) %>%
  mutate(chrom=c("1","2","3","4","5",
                 "6","7","8","9","10",
                 "11","12","13","14","15",
                 "16","17"))

#get threshold for significance in pvalue scale
-log10(max(cis_sQTL_sig$p.value))

sQTL_manhattan_plot <- ggplot(cis_sQTL,
                              aes(x=cum_pos, y=-log10(p.value),
                                  color=as.factor(chrom))) +
  
  ## mark inversion regions
  #geom_rect(data=inversion_regions,
  #          aes(xmin=inv_cumstart, xmax=inv_cumend,
  #              ymin=0,
  #              ymax=1.7),
  #          fill="red", inherit.aes = F) +
  
  # add each SNP–gene correlation.
  geom_point(shape=1) +
  scale_color_manual(values = rep(c("grey50", "black"),
                                  unique(length(axisdf_sQTL$chrom)))) +
  #scale_alpha_manual(name="FDR < .05",values = c(.1,1), labels=c("no", "yes")) +
  #scale_size_manual(name="FDR < .05", values = c(1,4), labels=c("no", "yes")) +
  #scale_color_gradientn(colors=c("grey","turquoise","blue","navyblue","black")) +
  
  geom_hline(lty=2, color="red", yintercept=4) +
  
  # custom axes
  scale_x_continuous(label = axisdf_sQTL$chrom,
                     breaks = axisdf_sQTL$center, expand=c(.01,.01)) +
  scale_y_continuous(expand=c(.01,.01)) +
  ylim(1.7,max(-log10(cis_sQTL$p.value))) +
  coord_cartesian(clip = 'off') +
  
  # custom theme
  theme_bw() +
  theme(legend.position = "none",
        text = element_text(size=24),
        axis.line.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x=element_text(size=18),
        plot.title = element_text(size=18),
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  labs(x="Chromosome", y="-log10 p-value", title = "cis-sQTL")

sQTL_manhattan_plot
eQTL_manhattan_plot

FIG_QTL <- eQTL_manhattan_plot / sQTL_manhattan_plot
ggsave("figures/FIG_QTL.png", plot = FIG_QTL, device = "png",
       width = 12, height=7.2, units = "in", dpi = 300)
