library(tidyr)
#library(gamlss) #for zero inflated beta regression
#library(Rfast2) #for zero inflated gamma regression?
library(glmmTMB)
library(performance)
library(DHARMa)
library(see)
library(dplyr)
library(tibble)

#### make a manhattan plot of genome-wide Fst ####
#Fst <- read.table("analysis/Fst/Fst_1Kb_windows.weir.txt", header=T) %>%
#  na.omit() %>%
#  rename(CHROM="chrom")
#Fst <- read.table("analysis/Fst/Fst_10Kb_windows.weir.txt", header=T) %>%
#  na.omit() %>%
#  rename("chrom"=CHROM)
Fst <- read.table("analysis/Fst/Fst_500kb_windows.weir.txt", header=T) %>%
  na.omit() %>%
  dplyr::rename(chrom=CHROM)

# convert negative Fst values to 0.
Fst$WEIGHTED_FST[Fst$WEIGHTED_FST < 0] <- 0
Fst$MEAN_FST[Fst$MEAN_FST < 0] <- 0

save(Fst, file = "data2/Rdata/Fst_500kb_windows.Rdata")

# get the cumulative length of each chromosome
chrom_lengths <- read.table("data/ref_genome_Ha412HO/chrom_sizes_Ha412HOv2.0.txt",
                            col.names = c("chrom", "length")) %>%
  mutate(length=as.numeric(length))

cum_lengths <- chrom_lengths %>%
  mutate(cumstart=cumsum(length)-length,
         cumend=cumsum(length))

# get inversion regions
inversion_regions_all <- read.table("analysis/inversions/inversion_regions.txt") %>%
  dplyr::select(V2,V3,V4,V9) %>%
  dplyr::rename(chrom=V2, start=V3, end=V4, name=V9) %>%
  #filter(name %in% c("pet05.01", "pet09.01", "pet11.01", "pet17.01")) %>%
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

Fst_manhattan_plot_supp <- ggplot(data=Fst_manhattan_df,
                             aes(x=cumstart, y=WEIGHTED_FST)) +
  # alternate shading of chromosomes
  geom_rect(data=cum_lengths, aes(x=NULL,y=NULL,
                                  xmin=cumstart, xmax=cumend,
                                  ymin=0,ymax=1,
                                  fill=as.factor(chrom)), alpha=0.4) +
  scale_fill_manual(values = rep(c("light grey", "white"), 17 )) +

  # add mean Fst line and inversion regions
  geom_line(linewidth=.125) +
  geom_rect(data=inversion_regions_all, aes(x=NULL, y=NULL,
                                        xmin=inv_cumstart, xmax=inv_cumend,
                                        ymin=0, ymax=1),
            fill="red",alpha=0.25) +
  # customize axes
  scale_x_continuous(label = axisdf_Fst$chrom, breaks = axisdf_Fst$center,
                     expand = c(.01,.01)) +
  scale_y_continuous(expand=c(.01,.01)) +
  coord_cartesian(clip = 'off') +
  
  # custom theme
  theme_bw(base_size=12) +
  theme(legend.position = "none",
        #axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),
        #text = element_text(size=24),
        #axis.text.x=element_blank(),
        #axis.line.x = element_blank(),
        axis.ticks = element_blank(),
        #axis.text.y = element_text(size=18),
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  labs(x="Chromosome", y=expression(italic(F)["st"]))
  
Fst_manhattan_plot_supp

ggsave("figures/Fst_manhattan_supp.pdf", device="pdf",plot=Fst_manhattan_plot_supp,
       width=175, height = 50, dpi = 300, units = "mm")

#### Regressions: LFC and dPSI vs Fst per gene ####

# get gene coordinates
expressed_genes <- read.table("analysis/DESeq2/expressed_genes.txt",
                              col.names = "Ha412_gene")
genes_gff <- read.table("data/ref_genome_Ha412HO/HAN412_Eugene_curated_v1_1.genes_gff.tmp") %>%
  dplyr::select(c(1:5,11)) #%>%
  #filter(!grepl("Chr00", V1)) #exclude ~550 genes on unplaced contigs

names(genes_gff) <- c("chrom","gene_start", "gene_end",
                      "gene_length", "strand", "Ha412_gene")

expressed_genes_gff <- genes_gff %>% filter(Ha412_gene %in% expressed_genes$Ha412_gene)

## set windows as the coordinates encompassing five sequential genes
#chrom <- ""
#for (i in 1:(length(genes_gff$Ha412_gene)-2)) {
#  if (genes_gff$chrom[i] == chrom) { #if gene is on same chrom as previous gene
#    if ( length(genes_gff$chrom[i-2])==0 ) { #if gene is the second gene of a new chrom
#      win_start <- 1 #genes_gff$gene_start[i-1]
#    } else if (genes_gff$chrom[i]!=genes_gff$chrom[i-2]) {
#      win_start <- 1 #genes_gff$gene_start[i-1]
#    } else { #it's the third or subsequent gene of the chrom
#      win_start <- genes_gff$gene_start[i-2]
#    }
#    if (genes_gff$chrom[i+2] != chrom) { #it's the second-to-last or last gene
#      if (genes_gff$chrom[i+1] != chrom) {#it's the last gene on chrom
#        win_end <- genes_gff$gene_end[i]
#      } else { #it's second to last
#        win_end <- genes_gff$gene_end[i+1]
#      }
#    } else { #it's not the last or second to last gene on the chrom
#      win_end <- genes_gff$gene_end[i+2]
#    }
#    
#  } else { #it's the first gene of a new chrom
#    chrom <- genes_gff$chrom[i]
#    win_start <- 1 #genes_gff$gene_start[i]
#    win_end <- genes_gff$gene_end[i+2]
#  }
# 
#  genes_gff$win_start[i] <- win_start
#  genes_gff$win_end[i] <- win_end
#}
#
## need to fill in the end positions of the last two genes
## bc our loop doesn't reach them
#genes_gff$win_start[nrow(genes_gff)-1] <- genes_gff$gene_start[nrow(genes_gff)-3]
#genes_gff$win_end[nrow(genes_gff)-1] <- genes_gff$gene_end[nrow(genes_gff)]
#genes_gff$win_start[nrow(genes_gff)] <- genes_gff$gene_start[nrow(genes_gff)-2]
#genes_gff$win_end[nrow(genes_gff)] <- genes_gff$gene_end[nrow(genes_gff)]
#tail(genes_gff)
#head(genes_gff)
#
#genes_gff <- genes_gff %>%
#  mutate(win_size = win_end - win_start + 1)
#
#write.table(dplyr::select(genes_gff, chrom, win_start, win_end),
#            file="analysis/Fst/custom_gene_windows.bed",
#            row.names = F, col.names = T, quote = F, sep = '\t')
#
## We will use this bed file to calculate Fst for the custom windows
## using Python scikit-allel (separate script)

# read in Fst gene windows results from python 
Fst_single_gene_windows <- read.csv("analysis/Fst/Fst_single_gene_windows.allel.csv",
                                    header=T) %>%
  dplyr::rename(gene_start=win_start, gene_end=win_end) %>%
  mutate(gene_start=gene_start+1) %>%
  left_join(genes_gff)
# convert negative Fst values to 0,
Fst_single_gene_windows$Fst[Fst_single_gene_windows$Fst < 0] <- 0

Fst_single_gene_windows_5kb_buffer <- read.csv("analysis/Fst/Fst_single_gene_windows_5kb_buffer.allel.csv",
                                    header=T) %>%
  mutate(win_start = win_start + 1,
         gene_start = win_start + 5000, gene_end = win_end - 5000) %>%
  left_join(genes_gff)
# convert negative Fst values to 0.
Fst_single_gene_windows_5kb_buffer$Fst[Fst_single_gene_windows_5kb_buffer$Fst < 0] <- 0

# load LFC and dPSI results
load("data2/Rdata/DESeq2_results_Shrink.Rdata") #from analyze_differential_expression.R
load("data2/Rdata/rmats_results_dfs.Rdata") #from analyze_splicing_rMATS.R


LFC_vs_Fst_df <- data.frame(de_results_Shrink) %>% 
  #filter(baseMean>=10) %>% 
  rownames_to_column("Ha412_gene") %>% 
  mutate(Ha412_gene=gsub(".*RNA","gene",Ha412_gene)) %>%
  left_join(Fst_single_gene_windows_5kb_buffer) %>% 
  dplyr::filter(!grepl("Chr00", Ha412_gene), grepl("gene:", Ha412_gene)) %>%
  na.omit() #only keep genes with values for both dPSI and Fst
# set LFC to zero if diff expression not significant
LFC_vs_Fst_df$log2FoldChange[LFC_vs_Fst_df$padj >= 0.05] <- 0

# <- all_AS_events_PSI %>%
#  na.omit() %>% 
#  dplyr::select(ID)
dPSI_vs_Fst_df <- all_AS_events_deltaPSI %>%
  #filter(ID %in% AS_events_no_miss$ID) %>%
  dplyr::rename(Ha412_gene=GeneID,deltaPSI=IncLevelDifference) %>%
  group_by(Ha412_gene) %>%
  slice_max(abs(deltaPSI)) %>%
  left_join(Fst_single_gene_windows_5kb_buffer) %>%
  dplyr::filter(!grepl("Chr00", Ha412_gene)) %>%
  as.data.frame(.) %>%
  na.omit() #only keep genes with values for both dPSI and Fst
# set dPSI to zero if not significant
dPSI_vs_Fst_df$deltaPSI[dPSI_vs_Fst_df$FDR >= 0.05] <- 0

# linear models
# zero-adjusted gamma regression ('ZAGA') for LFC
#fit_LFC_vs_Fst <- gamlss(abs(log2FoldChange) ~ Fst, family = ZAGA, data=LFC_vs_Fst_df)
#gamlss::Rsq(fit_LFC_vs_Fst, type = "both")
#plot(fit_LFC_vs_Fst)
fit_LFC_vs_Fst <- glmmTMB(abs(log2FoldChange) ~ Fst, family=ziGamma(link="log"), ziformula = ~1, data = LFC_vs_Fst_df)
tmp2 <- simulateResiduals(fit_LFC_vs_Fst)
performance::r2(fit_LFC_vs_Fst)
summary(fit_LFC_vs_Fst)
plot(tmp2)

# zero-inflated beta regression to model dPSI ~ Fst
#fit_dPSI_vs_Fst <- gamlss(abs(deltaPSI) ~ Fst, family = BEZI, data = dPSI_vs_Fst_df)
#plot(fit_dPSI_vs_Fst) #model diagnostics. looks okay.
#gamlss::Rsq(fit_dPSI_vs_Fst, type = "both")
fit_dPSI_vs_Fst <- glmmTMB(abs(deltaPSI) ~ Fst, family = beta_family(), ziformula = ~1, data=dPSI_vs_Fst_df)
tmp3 <- simulateResiduals(fit_dPSI_vs_Fst)
performance::r2(fit_dPSI_vs_Fst)
summary(fit_dPSI_vs_Fst)
plot(tmp3)

# get predictions from the models to plot trend lines
new_data_LFC <- data.frame(Fst=seq(0,1,.001))
#new_data_LFC$preds_LFC <- gamlss::predictAll(fit_LFC_vs_Fst, newdata = new_data_LFC,
#                                             type = "response")[[1]]
new_data_LFC$preds_LFC <- predict(object = fit_LFC_vs_Fst, newdata = new_data_LFC,
                                  type="response" )

new_data_dPSI <- data.frame(Fst=seq(0,1,.001))
#new_data_dPSI$preds_dPSI <- gamlss::predictAll(fit_dPSI_vs_Fst, newdata = new_data_dPSI,
#                                               type = "response")[[1]]
new_data_dPSI$preds_dPSI <- predict(object = fit_dPSI_vs_Fst, newdata = new_data_dPSI,
                                  type="response" )

# save dfs for plotting
save(LFC_vs_Fst_df, dPSI_vs_Fst_df,
     new_data_LFC, new_data_dPSI, file = "data2/Rdata/LFC_dPSI_vs_Fst.Rdata")

plot_LFC_vs_Fst <- ggplot(data=LFC_vs_Fst_df, aes(x=Fst, y=abs(log2FoldChange))) +
  geom_point(size=2.5,alpha=0.5) +
  geom_line(data=new_data_LFC, aes(x=Fst, y=preds_LFC),
            color="red", linewidth=1.5) +
  theme_bw(base_size = 8) +
  #ylim(c(0,10)) +
  theme(#text = element_text(size=24),#,
        #legend.title=element_blank(),
        #plot.title = element_text(size=18),
        #axis.text = element_text(size=18),
        panel.grid.major = element_blank()) +
  labs(x=expression(italic(F)["st"]), y="|log2FC|")

plot_dPSI_vs_Fst <- ggplot(data=dPSI_vs_Fst_df, aes(x=Fst, y=abs(deltaPSI))) +
  geom_point(size=2.5,alpha=0.5) +
  geom_line(data=new_data_dPSI, aes(x=Fst, y=preds_dPSI),
            color="red", linewidth=1.5) +
  theme_bw(base_size=8) +
  theme(#text = element_text(size=24),#,
        #legend.title=element_blank(),
      #plot.title = element_text(size=18),
      #axis.text = element_text(size=18),
        panel.grid.major = element_blank()) +
  labs(x=expression(italic(F)["st"]), y=~paste("|",Delta,"PSI|"))

#### Fst of DE/DS genes vs Fst of all genes ####
DE_DS_intersect_genes <- read.table("analysis/GO_analysis/study_intersect_DE–DS_genes.txt") %>%
  dplyr::rename(Ha412_gene=V1)
DE_only_genes <- read.table("analysis/GO_analysis/study_DE_genes_noLFCthreshold.txt") %>%
  dplyr::rename(Ha412_gene=V1) %>%
  filter(!Ha412_gene %in% DE_DS_intersect_genes$Ha412_gene)
DS_only_genes <- read.table("analysis/GO_analysis/study_DS_rMATS_genes.txt") %>%
  dplyr::rename(Ha412_gene=V1) %>%
  filter(!Ha412_gene %in% DE_DS_intersect_genes$Ha412_gene)
DE_DS_union_genes <- rbind(DE_DS_intersect_genes,
             rbind(DE_only_genes, DS_only_genes))
non_DE_DS_genes <- expressed_genes_gff %>%
  filter(!Ha412_gene %in% DE_DS_union_genes$Ha412_gene) %>%
  dplyr::select(Ha412_gene)
                                              

Fst_DE_only_genes <- Fst_single_gene_windows_5kb_buffer %>%
  filter(Ha412_gene %in% DE_only_genes$Ha412_gene) %>%
  na.omit() %>%
  mutate(type="DE_only")
Fst_DS_only_genes <- Fst_single_gene_windows_5kb_buffer %>%
  filter(Ha412_gene %in% DS_only_genes$Ha412_gene) %>%
  na.omit() %>%
  mutate(type="DS_only")
Fst_DE_DS_intersect_genes <- Fst_single_gene_windows_5kb_buffer %>%
  filter(Ha412_gene %in% DE_DS_intersect_genes$Ha412_gene) %>%
  na.omit() %>%
  mutate(type="DE_DS_intersect")
Fst_non_DE_DS_genes <- Fst_single_gene_windows_5kb_buffer %>%
  filter(Ha412_gene %in% non_DE_DS_genes$Ha412_gene) %>%
  na.omit() %>%
  mutate(type="non_DE_DS")

mean(Fst_DE_only_genes$Fst)
mean(Fst_DS_only_genes$Fst)
mean(Fst_DE_DS_intersect_genes$Fst)
mean(Fst_non_DE_DS_genes$Fst)
#mean(Fst_single_gene_windows_5kb_buffer$Fst, na.rm=T)
Fst_genes_df <- full_join(Fst_DE_only_genes, Fst_DS_only_genes) %>%
  full_join(Fst_DE_DS_intersect_genes) %>%
  full_join(Fst_non_DE_DS_genes) %>%
  dplyr::select(Ha412_gene, Fst, type)
Fst_genes_df$type <- factor(Fst_genes_df$type, levels=c( "non_DE_DS","DE_only",
                                                         "DS_only",
                                                         "DE_DS_intersect"))

kruskal.test(Fst ~ type, data = Fst_genes_df)
pairwise.wilcox.test(Fst_genes_df$Fst, Fst_genes_df$type,
                     p.adjust.method = "BH")

save(Fst_genes_df, file = "data2/Rdata/Fst_genes_df.Rdata")
plot_Fst_vs_gene_sets <- ggplot(data=Fst_genes_df, aes(x=type, y=Fst, fill=type)) +
  geom_boxplot(alpha=0.5,linewidth=.25, outlier.size = .5, outlier.alpha = .25) +
  scale_fill_manual(values = c("grey","#317EC2","red","purple")) +
  scale_x_discrete(labels=c("Non-\nDE/DS", "DE", "DS","DE∩DS")) +
  scale_y_continuous(limits = c(0,1.05)) +
  theme_bw(base_size=8) +
  theme(legend.position = "none") +
  labs(x="Gene set",y=expression(italic(F)["st"]))
ggsave("figures/Fst_gene_set_boxplot.pdf", plot=last_plot(), device = "pdf",
       width = 58.33, height = 43.75, dpi = 300)

#### Fst of splice site variants ####
library(GenomicRanges)
splice_var_Fst <- read.table("analysis/snpEff/splice_variants_Fst.txt") %>%
  dplyr::rename(chrom=V1,pos=V2,Fst=V3)
# convert negative Fst values to 0,
splice_var_Fst$Fst[splice_var_Fst$Fst < 0] <- 0

mean(splice_var_Fst$Fst) #mean Fst = 0.087

# compare to all Fst of all SNPs
per_site_Fst <- read.table("analysis/Fst/Fst_per_site.weir.txt", header=T) %>%
  na.omit()
# convert negative Fst values to 0.
per_site_Fst$WEIR_AND_COCKERHAM_FST[per_site_Fst$WEIR_AND_COCKERHAM_FST < 0] <- 0
mean(per_site_Fst$WEIR_AND_COCKERHAM_FST) #mean transcriptome-wide Fst = 0.084

# save list of splice vars with highest Fst
top_splice_var_Fst <- splice_var_Fst %>%
  #filter(Fst > mean(per_site_Fst$WEIR_AND_COCKERHAM_FST) +
  #         4*sd(per_site_Fst$WEIR_AND_COCKERHAM_FST)) %>%
  filter(Fst >  quantile(per_site_Fst$WEIR_AND_COCKERHAM_FST, .95)) %>%
  mutate(pos_0based = pos - 1) %>%
  dplyr::select(chrom, pos_0based, pos, Fst)
write.table(top_splice_var_Fst, "analysis/snpEff/top_splice_variants_Fst.txt",
            row.names = F, quote = F, sep = "\t")

#### Fst of spliceosomal genes ####
spliceosomal_genes_prot <- 
  read.table("analysis/BLAST_out/spliceosome_components_prot_vs_HAN412_1e-20.txt") %>%
  dplyr::rename(blast_hit=V1,Ha412_gene=V2) %>%
  mutate(Ha412_gene=gsub("mRNA","gene",Ha412_gene)) %>%
  dplyr::select(Ha412_gene) %>%
  unique() %>%
  left_join(dplyr::select(Fst_single_gene_windows_5kb_buffer, Ha412_gene, Fst)) %>%
  na.omit()

spliceosomal_genes_ncRNA <- 
  read.table("analysis/BLAST_out/spliceosome_components_ncRNA_vs_HAN412_1e-20.txt") %>%
  dplyr::rename(blast_hit=V1,Ha412_gene=V2) %>%
  mutate(Ha412_gene=gsub("mRNA","gene",Ha412_gene)) %>%
  dplyr::select(Ha412_gene) %>%
  unique() %>%
  left_join(dplyr::select(Fst_single_gene_windows_5kb_buffer, Ha412_gene, Fst)) %>%
  na.omit()

spliceosomal_genes <- rbind(spliceosomal_genes_ncRNA, spliceosomal_genes_prot) %>%
  unique()
top_spliceosomal_genes_Fst <- spliceosomal_genes %>% 
  filter(Fst > quantile(Fst_single_gene_windows_5kb_buffer$Fst, .95, na.rm=T)) %>%
  left_join(Ha412_Ath_mappings)
  #filter(Fst > mean(Fst_single_gene_windows$Fst, na.rm = T) + 
  #         4*sd(Fst_single_gene_windows$Fst, na.rm = T))
write.table(top_spliceosomal_genes_Fst, "analysis/Fst/top_spliceosomal_genes_Fst.txt",
            row.names = F, quote = F, sep = "\t")

mean(spliceosomal_genes$Fst)
mean(Fst_single_gene_windows_5kb_buffer$Fst, na.rm = T)
