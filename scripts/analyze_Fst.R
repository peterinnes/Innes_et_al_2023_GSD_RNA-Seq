library(tidyr)

#### make a manhattan plot of genome-wide Fst ####
#Fst <- read.table("analysis/Fst/Fst_1Kb_windows.weir.txt", header=T) %>%
#  na.omit() %>%
#  rename(CHROM="chrom")
#Fst <- read.table("analysis/Fst/Fst_10Kb_windows.weir.txt", header=T) %>%
#  na.omit() %>%
#  rename("chrom"=CHROM)
Fst <- read.table("analysis/Fst/Fst_2Mb_windows.weir.txt", header=T) %>%
  na.omit() %>%
  rename(chrom=CHROM)

# convert negative Fst values to 0.
for (i in 1:length(Fst$WEIGHTED_FST)) {
  if (Fst$WEIGHTED_FST[i] < 0) {
   Fst$WEIGHTED_FST[i] <- 0
  }
  if (Fst$MEAN_FST[i] < 0) {
     Fst$MEAN_FST[i] <- 0
  }
}

save(Fst, file = "data2/Rdata/Fst_2Mb_windows.Rdata")

# get the cumulative position of each SNP/window and then the max SNP end position on each chrom
# get the cumulative length of each chromosome
chrom_lengths <- read.table("data/ref_genome_Ha412HO/chrom_sizes_Ha412HOv2.0.txt",
                            col.names = c("chrom", "length")) %>%
  mutate(length=as.numeric(length))

cum_lengths <- chrom_lengths %>%
  mutate(cumstart=cumsum(length)-length,
         cumend=cumsum(length))

# get inversion regions
inversion_regions <- read.table("analysis/inversions/inversion_regions.txt") %>%
  dplyr::select(V2,V3,V4,V9) %>%
  rename(chrom=V2, start=V3, end=V4, name=V9) %>%
  filter(name %in% c("pet05.01", "pet09.01", "pet11.01", "pet17.01")) %>%
  left_join(cum_lengths[c(1,3)]) %>%
  mutate(inv_cumstart=start+cumstart, inv_cumend=end+cumstart)

Fst_manhattan_df <- Fst %>% 
  inner_join(cum_lengths, by = "chrom") %>% 
  mutate(cumstart = BIN_START + cumstart) %>%
  filter(!grepl("Chr00", chrom))
#Fst$CHROM <- as.factor(Fst$CHROM)

# get mid point of each chromosome
axisdf_Fst <- Fst_manhattan_df %>%
  group_by(chrom) %>%
  summarize(center=(max(cumstart) + min(cumstart)) / 2) %>%
  mutate(chrom=c("1","2","3","4","5",
                 "6","7","8","9","10",
                 "11","12","13","14","15",
                 "16","17"))

Fst_manhattan_plot <- ggplot(data=Fst_manhattan_df,
                             aes(x=cumstart, y=WEIGHTED_FST)) +
  # alternate shading of chromosomes
  geom_rect(data=cum_lengths, aes(x=NULL,y=NULL,
                                  xmin=cumstart, xmax=cumend,
                                  ymin=0,ymax=.8,
                                  fill=as.factor(chrom)), alpha=0.4) +
  scale_fill_manual(values = rep(c("light grey", "white"), 17 )) +

  # add mean Fst line and inversion regions
  geom_line() +
  geom_rect(data=inversion_regions, aes(x=NULL, y=NULL,
                                        xmin=inv_cumstart, xmax=inv_cumend,
                                        ymin=0, ymax=0.8),
            fill="red",alpha=0.25) +
  # customize axes
  scale_x_continuous(label = axisdf_Fst$chrom, breaks = axisdf_Fst$center,
                     expand = c(.01,.01)) +
  scale_y_continuous(expand=c(.01,.01)) +
  coord_cartesian(clip = 'off') +
  
  # custom theme
  theme_bw(base_size=18) +
  theme(legend.position = "none",
        #axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),
        #text = element_text(size=24),
        axis.text.x=element_blank(),
        axis.line.x = element_blank(),
        axis.ticks = element_blank(),
        #axis.text.y = element_text(size=18),
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  labs(x="", y=expression(italic(F)["st"]))
  
Fst_manhattan_plot

ggsave("figures/Fst_manhattan.png", device="png",plot=Fst_manhattan_plot,
       width=8, height = 6, dpi = 300, units = "in")

#### expr and splicing vs Fst per gene (windows around gene) ####

# make custom windows around genes
expressed_genes <- read.table("analysis/DESeq2/expressed_genes.txt",
                              col.names = "Ha412_gene")
genes_gff <- read.table("data/ref_genome_Ha412HO/HAN412_Eugene_curated_v1_1.genes_gff.tmp") %>%
  dplyr::select(c(1:5,11)) %>%
  filter(!grepl("Chr00", V1)) #exclude ~550 genes on unplaced contigs

names(genes_gff) <- c("chrom","gene_start", "gene_end",
                      "gene_length", "strand", "Ha412_gene")

genes_gff <- genes_gff %>% filter(Ha412_gene %in% expressed_genes$Ha412_gene)

# set windows as the coordinates encompassing five sequential genes
chrom <- ""
for (i in 1:(length(genes_gff$Ha412_gene)-2)) {
  if (genes_gff$chrom[i] == chrom) { #if gene is on same chrom as previous gene
    if ( length(genes_gff$chrom[i-2])==0 ) { #if gene is the second gene of a new chrom
      win_start <- 1 #genes_gff$gene_start[i-1]
    } else if (genes_gff$chrom[i]!=genes_gff$chrom[i-2]) {
      win_start <- 1 #genes_gff$gene_start[i-1]
    } else { #it's the third or subsequent gene of the chrom
      win_start <- genes_gff$gene_start[i-2]
    }
    if (genes_gff$chrom[i+2] != chrom) { #it's the second-to-last or last gene
      if (genes_gff$chrom[i+1] != chrom) {#it's the last gene on chrom
        win_end <- genes_gff$gene_end[i]
      } else { #it's second to last
        win_end <- genes_gff$gene_end[i+1]
      }
    } else { #it's not the last or second to last gene on the chrom
      win_end <- genes_gff$gene_end[i+2]
    }
    
  } else { #it's the first gene of a new chrom
    chrom <- genes_gff$chrom[i]
    win_start <- 1 #genes_gff$gene_start[i]
    win_end <- genes_gff$gene_end[i+2]
  }
 
  genes_gff$win_start[i] <- win_start
  genes_gff$win_end[i] <- win_end
}

# need to fill in the end positions of the last two genes
# bc our loop doesn't reach them
genes_gff$win_start[nrow(genes_gff)-1] <- genes_gff$gene_start[nrow(genes_gff)-3]
genes_gff$win_end[nrow(genes_gff)-1] <- genes_gff$gene_end[nrow(genes_gff)]
genes_gff$win_start[nrow(genes_gff)] <- genes_gff$gene_start[nrow(genes_gff)-2]
genes_gff$win_end[nrow(genes_gff)] <- genes_gff$gene_end[nrow(genes_gff)]
tail(genes_gff)
head(genes_gff)

genes_gff <- genes_gff %>%
  mutate(win_size = win_end - win_start + 1)

write.table(dplyr::select(genes_gff, chrom, win_start, win_end),
            file="analysis/Fst/custom_gene_windows.bed",
            row.names = F, col.names = T, quote = F, sep = '\t')

# We will use this bed file to calculate Fst for the custom windows
# using Python scikit-allel (separate script)

# read in Fst custom gene windows results from python 
Fst_custom_gene_windows <- read.csv("analysis/Fst/Fst_custom_gene_windows.allel.csv",
                                    header=T) %>%
  left_join(genes_gff)

# convert negative Fst values to 0,
# and missing values (windows with no variants) to 0.
for (i in 1:nrow(Fst_custom_gene_windows) ){
  if ( is.na(Fst_custom_gene_windows$Fst[i]) ){
    Fst_custom_gene_windows$Fst[i] <- 0
  } else if (Fst_custom_gene_windows$Fst[i] < 0) {
    Fst_custom_gene_windows$Fst[i] <- 0
  }
}
  
head(Fst_custom_gene_windows)

# load LFC and dPSI results
load("data2/Rdata/DESeq2_results_Shrink.Rdata") #from analyze_differential_expression.R
load("data2/Rdata/rmats_results_dfs.Rdata") #from analyze_splicing_rMATS.R

LFC_vs_Fst_df <- data.frame(de_results_Shrink) %>% # without stringtie
  rownames_to_column("Ha412_gene") %>% 
  mutate(Ha412_gene=gsub(".*RNA","gene",Ha412_gene)) %>%
  left_join(Fst_custom_gene_windows) %>% 
  dplyr::filter(!grepl("Chr00", Ha412_gene), grepl("gene:", Ha412_gene))
# set LFC to zero if not significant 
for (i in 1:nrow(LFC_vs_Fst_df) ){
  if ( !is.na(LFC_vs_Fst_df$padj[i])){
    if ( LFC_vs_Fst_df$padj[i] >= .05 ){
      LFC_vs_Fst_df$log2FoldChange[i] <- 0
    }
  }
}


# with splicing data
dPSI_vs_Fst_df <- all_AS_events_deltaPSI
# set dPSI to zero if not significant
for (i in 1:nrow(dPSI_vs_Fst_df) ){
  if (dPSI_vs_Fst_df$FDR[i] >= .05 ){
    dPSI_vs_Fst_df$IncLevelDifference[i] <- 0
  }
}
dPSI_vs_Fst_df <- dPSI_vs_Fst_df%>%
  rename(Ha412_gene=GeneID) %>%
  group_by(Ha412_gene) %>%
  summarise(n_events=n(),max_deltaPSI=max(abs(IncLevelDifference))) %>%
  left_join(Fst_custom_gene_windows) %>%
  dplyr::filter(!grepl("Chr00", Ha412_gene)) %>%
  as.data.frame()

# save dfs for plotting
save(LFC_vs_Fst_df, dPSI_vs_Fst_df, file = "data2/Rdata/LFC_dPSI_vs_Fst.Rdata")

# plotting and linear models
plot_LFC_vs_Fst <- ggplot(data=LFC_vs_Fst_df, aes(x=Fst, y=abs(log2FoldChange))) +
  geom_point(size=2.5,alpha=0.5) +
  geom_smooth(method = "lm", color="red") +
  theme_bw(base_size = 18) +
  #ylim(c(0,10)) +
  theme(#text = element_text(size=24),#,
        #legend.title=element_blank(),
        #plot.title = element_text(size=18),
        #axis.text = element_text(size=18),
        panel.grid.major = element_blank()) +
  labs(x=expression(italic(F)["st"]), y="|log2FC|")

fit_LFC_vs_Fst <- lm(abs(log2FoldChange) ~ Fst, data = LFC_vs_Fst_df)
summary(fit_LFC_vs_Fst)

plot_dPSI_vs_Fst <- ggplot(data=dPSI_vs_Fst_df, aes(x=Fst, y=deltaPSI)) +
  geom_point(size=2.5,alpha=0.5) +
  geom_smooth(method = "lm", color="red") +
  theme_bw(base_size=18) +
  theme(#text = element_text(size=24),#,
        #legend.title=element_blank(),
      #plot.title = element_text(size=18),
      #axis.text = element_text(size=18),
        panel.grid.major = element_blank()) +
  labs(x=expression(italic(F)["st"]), y=~paste("|",Delta,"PSI|"))

fit_dPSI_vs_Fst <- lm(max_deltaPSI ~ Fst, data = dPSI_vs_Fst_df)
summary(fit_dPSI_vs_Fst)


#### scraps ####
## Using windowed Fst from vcftools
#Fst_per_gene_window_df <- read.table("analysis/Fst/intersect_Fst_100Kb_windows_gff_genes.txt") %>% 
#  dplyr::select(c(1:6, 8:10, 17))
#names(Fst_per_gene_window_df) <- c("CHROM","BIN_START","BIN_END","N_VARIANTS",
#                                   "WEIGHTED_FST", "MEAN_FST","GENE_START",
#                                   "GENE_END", "GENE_LENGTH", "Ha412_gene")
#for (i in 1:length(Fst_per_gene_window_df$WEIGHTED_FST)) {
#  if (Fst_per_gene_window_df$WEIGHTED_FST[i] < 0) {
#    Fst_per_gene_window_df$WEIGHTED_FST[i] <- 0
#  }
#  if (Fst_per_gene_window_df$MEAN_FST[i] < 0) {
#    Fst_per_gene_window_df$MEAN_FST[i] <- 0
#  }
#}
#
#tmp <- Fst_per_gene_window_df %>% group_by(Ha412_gene) %>%
#  summarise("Fst"=mean(WEIGHTED_FST))

## make windows of certain bp size around each gene
#genes_gff <- genes_gff %>%
#  mutate(window_start=floor( (gene_start+(gene_length/2)) - 250000 ),
#         window_end=ceiling( (gene_start+(gene_length/2)) + 250000 ))
#
#for (i in 1:length(genes_gff$Ha412_gene)) {
#  if (genes_gff$window_start[i] < 1) {
#    genes_gff$window_start[i] <- 1
#  }
#}
#




