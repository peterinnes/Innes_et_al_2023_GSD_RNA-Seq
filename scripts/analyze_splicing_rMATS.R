# analyzing rMATS output #
BiocManager::install("preprocessCore")
library(dplyr)
library(tidyr)
library(ggplot2)
library(preprocessCore)

#### check the results summary ####
# How many significant DS events of the 5 different types of splicing?
rmatsout_summary <- read.table('analysis/rMATS/results_2022-07-14/summary.txt', header = T)
rmatsout_summary <- pivot_longer(rmatsout_summary,2:9, names_to = "Category", values_to = "Count") %>%
  subset(Category %in% c("TotalEventsJCEC", "SignificantEventsJCEC"))

rmatsout_summary$EventType <- factor(rmatsout_summary$EventType,
                                     levels=c("RI", "A3SS", "A5SS", "MXE", "SE"))

# plot non-significant events and significant events together
plot_rmats_event_counts <- ggplot(data = rmatsout_summary %>% subset(Category=="TotalEventsJCEC"),
       aes(x=EventType, y=Count)) +
  geom_col(position = 'identity', color="grey50", alpha=0.5) +
  geom_col(data=rmatsout_summary %>% subset(Category=="SignificantEventsJCEC"),
           aes(x=EventType, y=Count),fill="black", position = 'identity') +
  theme_bw() +
  theme(text = element_text(size=24),
        axis.text = element_text(size=18),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  labs(x="Event Type")

ggsave("figures/plot_rMATS_event_counts_raw.png", plot=plot_rmats_event_counts,
       device = "png",
       height = 6, width = 6, dpi = 300, units = "in")

# significant events bar plot
ggplot(data = rmatsout_summary %>% subset(Category=="SignificantEventsJCEC"),
       aes(x=EventType, y=Count)) +
  geom_bar(stat = 'identity') +
  theme_classic() +
  theme(text = element_text(size=24),
        axis.text = element_text(size=18))
  
  

#### merge all the results files together ####
# add splice-type prefix to all the event IDs so we can tell different event types apart after we merge them all together
SE_events <- read.table('analysis/rMATS/results_2022-07-14/SE.MATS.JCEC.txt', header = T) %>%
  mutate(ID = gsub("^", "SE_", ID))
A5SS_events <- read.table('analysis/rMATS/results_2022-07-14/A5SS.MATS.JCEC.txt', header = T) %>%
  mutate(ID = gsub("^", "A5SS_", ID))
A3SS_events <- read.table('analysis/rMATS/results_2022-07-14/A3SS.MATS.JCEC.txt', header = T) %>% 
  mutate(ID = gsub("^", "A3SS_", ID))
MXE_events <- read.table('analysis/rMATS/results_2022-07-14/MXE.MATS.JCEC.txt', header = T) %>%
  mutate(ID = gsub("^", "MXE_", ID))
RI_events <- read.table('analysis/rMATS/results_2022-07-14/RI.MATS.JCEC.txt', header = T) %>% 
  mutate(ID = gsub("^", "RI_", ID))

# count number of significant intron retention splicing events, Should match the summary
RI_events |>
  subset(FDR < .05) |>
  #subset(abs(IncLevelDifference) > .05) |>
  dim() 

# combine all the results files into one file, for use in PCA of the Inclusion Levels, in the script pca.R
all_AS_events <- full_join(SE_events, RI_events) %>%
  full_join(MXE_events) %>%
  full_join(A5SS_events) %>% 
  full_join(A3SS_events) %>%
  arrange(FDR)

all_AS_events_PSI <- all_AS_events %>%
  dplyr::select(ID, GeneID, chr, IncLevel1, IncLevel2) %>%
  separate(IncLevel1, into = c("Kane-602-10_S8_R1_001.trimmed.fq.gz",
                               "Kane-602-11_S9_R1_001.trimmed.fq.gz",
                               "Kane-602-12_S10_R1_001.trimmed.fq.gz",
                               "Kane-602-14_S11_R1_001.trimmed.fq.gz",
                               "Kane-602-15_S12_R1_001.trimmed.fq.gz",
                               "Kane-602-1_S1_R1_001.trimmed.fq.gz",
                               "Kane-602-2_S2_R1_001.trimmed.fq.gz",
                               "Kane-602-3_S3_R1_001.trimmed.fq.gz",
                               "Kane-602-4_S4_R1_001.trimmed.fq.gz",
                               "Kane-602-7_S5_R1_001.trimmed.fq.gz",
                               "Kane-602-8_S6_R1_001.trimmed.fq.gz",
                               "Kane-602-9_S7_R1_001.trimmed.fq.gz"), sep = ",") %>%
  separate(IncLevel2, into = c("Kane-602-16_S13_R1_001.trimmed.fq.gz",
                               "Kane-602-17_S14_R1_001.trimmed.fq.gz",
                               "Kane-602-18_S15_R1_001.trimmed.fq.gz",
                               "Kane-602-19_S16_R1_001.trimmed.fq.gz",
                               "Kane-602-21_S17_R1_001.trimmed.fq.gz",
                               "Kane-602-23_S18_R1_001.trimmed.fq.gz",
                               "Kane-602-24_S19_R1_001.trimmed.fq.gz",
                               "Kane-602-25_S20_R1_001.trimmed.fq.gz",
                               "Kane-602-26_S21_R1_001.trimmed.fq.gz",
                               "Kane-602-27_S22_R1_001.trimmed.fq.gz",
                               "Kane-602-28_S23_R1_001.trimmed.fq.gz",
                               "Kane-602-30_S24_R1_001.trimmed.fq.gz"), sep = ",")

#get_log2FC <- function(A,B){
#  log2FC <- log2(A) - log2(B)
#  return(log2FC)
#}
## calculate log2FoldChange for dune non-dune 
#all_AS_events_PSI$log2FC <- get_log2FC(rowMeans(all_AS_events_PSI[4:15], na.rm = T),
#                                       rowMeans(all_AS_events_PSI[16:27], na.rm = T))

## calculate deltaPSI for dune - non-dune
#all_AS_events_PSI$deltaPSI <- rowMeans(all_AS_events_PSI[4:15],
#na.rm = T) - rowMeans(all_AS_events_PSI[16:27], na.rm = T)

write.table(all_AS_events_PSI, "analysis/rMATS/results_2022-07-14/all_AS_events_PSI.txt", row.names = F, quote = F, sep = '\t')

sig_AS_events <- subset(all_AS_events, FDR<.05) %>%
  rename("Ha412_gene"=GeneID)

#### make "phenotype" table for sQTL and (?) PCA?? a la leafcutter sQTL method ####
#rmats_splice_pheno_df <- all_AS_events_PSI[,4:27] %>%
#  lapply(as.numeric) %>% as.data.frame()
#rownames(rmats_splice_pheno_df) <- all_AS_events_PSI$ID
#
## throw out splice event if > 40% missingness
#tmp <- rmats_splice_pheno_df %>%
#  mutate(missing_perc = rowMeans(is.na(rmats_splice_pheno_df))) %>%
#  filter(missing_perc<.4) %>%
#  dplyr::select(!missing_perc)
#         
## set missing values to the average PSI of each event (each row)
#for( i in 1:nrow(tmp) ){
#  for( j in 1:ncol(tmp) ){
#    if( is.na(tmp[i,j]) ){
#      tmp[i,j] <- rowMeans(tmp[i,], na.rm = T)
#    }
#  }
#}
#
## if there is too little variation, remove
#tmp <- transform(tmp, sd=apply(tmp,1, sd, na.rm = TRUE)) %>%
#  filter(sd>=.005) %>%
#  dplyr::select(!sd)
#
## standardize across individuals
## have to transpose the df because scale() operates within columns
#tmp_scaled <- scale(t(tmp))
## transpose back again to quantile normalize across splice events
#rmats_splice_pheno_df <- data.frame(normalize.quantiles(t(tmp_scaled)))
#rownames(rmats_splice_pheno_df) <- rownames(tmp)
#colnames(rmats_splice_pheno_df) <- colnames(tmp) %>%
#  #gsub("_R1_001.trimmed.fq.gz","",.) %>%
#  gsub("Kane\\.602\\.","Kane-602-",.)
## reorder columns to be consistent with order in vcf file
#rmats_splice_pheno_df <- rmats_splice_pheno_df %>% 
#  dplyr::select(as.factor(sample_table_deseq2$fastq))


## filter out low-expressed-genes / genes with low support?
#low_expressed_genes$groupID
#all_AS_events <- all_AS_events %>%
#  mutate(groupID = gsub("gene:", "", GeneID))
#rmats_keep <- !(all_AS_events$groupID %in% low_expressed_genes$groupID) # check if any AS genes are in low-expressed list, #then flip the boolean to make a list of genes to keep. This only filters 110 genes, so the built-in rMATS filter must have #gotten most of them? This step probably isn't worth it. 
#temp <- all_AS_events[keep,]

## trying a different filter: at least 5 reads for an event in at least 3 different individuals, across inclusion/skip form regardless of habitat. The default rMATS filter keeps all events with at least one read supporting the exon inclusion form and the exon skipping form in either dune or non-dune samples. Put another way, The filter removes events where at least one of the sample groups does not have any reads to support the event (no inclusion reads and no skipping reads). The filter also removes events where neither sample group has a read for one of the isoforms (neither sample group has inclusion or neither sample group has skipping).

## first, split the columns with per event reads counts of each category (habitat*inclusion)
#temp <- all_AS_events %>%
#  `rownames<-`(all_AS_events$ID) %>%
#  dplyr::select(IJC_SAMPLE_1, IJC_SAMPLE_2, SJC_SAMPLE_1, SJC_SAMPLE_2) %>%
#  separate(IJC_SAMPLE_1, sep = ',', into = paste0("IJC1_",dune_fq)) %>% 
#  separate(IJC_SAMPLE_2, sep = ',', into = paste0("IJC2_",non.dune_fq)) %>% 
#  separate(SJC_SAMPLE_1, sep = ',', into = paste0("SJC1_",dune_fq)) %>% 
#  separate(SJC_SAMPLE_2, sep = ',', into = paste0("SJC2_",non.dune_fq)) %>%
#  mutate(across(everything(), as.numeric))
#
##temp_keep_IJC1 <- rowSums(temp[1:12] >= 5) >= 3
##temp_keep_IJC2 <- rowSums(temp[13:24] >= 5) >= 3
##temp_keep_SJC1 <- rowSums(temp[25:36] >= 5) >= 3
##temp_keep_SJC2 <- rowSums(temp[37:48] >= 5) >= 3

#temp_keep_IJC <- rowSums(temp[1:24] >= 5) >= 3 # for each splice event, require at least 3 individuals, irregardless of habitat, to have at least 5 reads mapping to the Inclusion form.  
#temp_keep_SJC <- rowSums(temp[25:48] >= 5) >= 3

## intersect the vectors for TRUE
##rmats_keep <- temp_keep_IJC1 & temp_keep_IJC2 & temp_keep_SJC1 & temp_keep_SJC2
#rmats_keep <- temp_keep_IJC & temp_keep_SJC
#sum(!rmats_keep) #this filters out 86276 events, leaving 9657
#all_AS_events_filtered <- (all_AS_events[rmats_keep,])

#### Splicing manhattan plot ####
all_AS_events_deltaPSI <- all_AS_events %>%
  mutate(exon_start = coalesce(exonStart_0base, X1stExonStart_0base,
                               longExonStart_0base, riExonStart_0base),
         exon_end = coalesce(exonEnd, X1stExonEnd,
                             longExonEnd, riExonEnd)) %>%
  dplyr::select(ID, GeneID, "chrom"=chr, exon_start, exon_end, IncLevelDifference, FDR) %>%
  mutate(chrom=gsub("chr", "", chrom)) %>%
  filter(ID %in% rownames(rmats_splice_pca_df)) #filter out events with > 40% missingness


#for( i in 1:nrow(all_AS_events_deltaPSI)){
#  if( all_AS_events_deltaPSI$FDR[i] == 0 ){
#    all_AS_events_deltaPSI$FDR[i] <- 1.134807e-14 #next smallest FDR
#  }
#}

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
  rename("chrom"=V2, "start"=V3, "end"=V4, "name"=V9) %>%
  filter(name %in% c("pet05.01", "pet09.01", "pet11.01", "pet17.01")) %>%
  left_join(cum_lengths[c(1,3)]) %>%
  mutate(inv_cumstart=start+cumstart, inv_cumend=end+cumstart)

splice_manhattan_df <- inner_join(all_AS_events_deltaPSI, cum_lengths) %>%
  #rowwise() %>%
  mutate(exon_startcum = exon_start + cumstart, sig=if_else(FDR<.05, "yes", "no"))

axisdf_splice <- splice_manhattan_df %>%
  group_by(chrom) %>%
  summarize(center=(max(exon_startcum) + min(exon_startcum)) / 2) %>%
  mutate(chrom=c("1","2","3","4","5",
                 "6","7","8","9","10",
                 "11","12","13","14","15",
                 "16","17"))

# split apart significant DE genes and nonsignifcant DE genes so we can downsample the nonsig data, for better vis
sig_splice_data <- splice_manhattan_df %>% 
  subset(FDR < 0.05)
notsig_splice_data <- splice_manhattan_df %>% 
  subset(FDR >= 0.05) %>%
  group_by(chrom) %>% 
  sample_frac(0.5)

splice_manhattan_df_reduced <- bind_rows(sig_splice_data, notsig_splice_data) 

splice_manhattan_plot <- ggplot(splice_manhattan_df_reduced,
                                aes(x=exon_startcum, y=IncLevelDifference),
                                alpha=IncLevelDifference) +
                                    #color=-log10(FDR))) +
  # alternate shading of chromosomes
  geom_rect(data=cum_lengths,
            aes(xmin=cumstart, xmax=cumend,
                ymin=-1,
                ymax=1,
                fill=as.factor(chrom)),
            inherit.aes = F,
            alpha=0.5) +
  scale_fill_manual(values = rep(c("light grey", "white"), 17), guide = "none") +
  
  # mark inversion regions
  geom_rect(data=inversion_regions,
            aes(xmin=inv_cumstart, xmax=inv_cumend,
                ymin=-1,
                ymax=1),
            fill="red", alpha=0.25, inherit.aes = F) +
  
  # add each splice event
  geom_point(aes(size=as.factor(sig=="no")), shape=1) +
  #scale_alpha_manual(name="FDR < .05",values = c(.1,.5), labels=c("no", "yes")) +
  scale_size_manual(name="Significance (FDR < .05)", values = c(4,1), labels=c("yes", "no")) +
  #geom_point(alpha=0.75, aes(size=-log(FDR))) +
  #scale_color_gradientn(colors=c("grey","turquoise","blue","navyblue","black")) +
  
  # custom axes
  scale_x_continuous(label = axisdf_splice$chrom,
                     breaks = axisdf_splice$center, expand=c(0.01,0.01)) +
  scale_y_continuous(expand=c(0.01,0.01)) +
  
  # custom theme
  theme_bw() +
  theme(legend.position = "top",
        #axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),
        text = element_text(size=24),
        axis.line.x = element_blank(),
        axis.ticks = element_blank(),
        #axis.text.x=element_blank(),
        #plot.title = element_text(size=18),
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.text = element_text(size=18),
        legend.title = element_text(size=18),
        axis.text = element_text(size=18),
        legend.box.margin = margin(t = -30)) +
  
  labs(x="Chromosome",y=~paste(Delta,"PSI"))
splice_manhattan_plot


## alt version with log2FC instead of deltaPSI
#
#splice_manhattan_plot_alt <- ggplot(splice_manhattan_df_reduced,
#                                aes(x=exon_startcum, y=log2FC,
#                                    color=-log10(FDR))) +
#  # alternate shading of chromosomes
#  geom_rect(data=cum_lengths,
#            aes(xmin=cumstart, xmax=cumend,
#                ymin=min(splice_manhattan_df$log2FC, na.rm = T),
#                ymax=max(splice_manhattan_df$log2FC, na.rm = T),
#                fill=as.factor(chrom)),
#            inherit.aes = F,
#            alpha=0.5) +
#  scale_fill_manual(values = rep(c("light grey", "white"), 17)) +
#  
#  # mark inversion regions
#  geom_rect(data=inversion_regions,
#            aes(xmin=inv_cumstart, xmax=inv_cumend,
#                ymin=min(splice_manhattan_df$log2FC, na.rm = T),
#                ymax=max(splice_manhattan_df$log2FC, na.rm = T)),
#            fill="red", alpha=0.25, inherit.aes = F) +
#  
#  # add each splice event
#  geom_point(alpha=0.75, aes(size=-log(FDR))) +
#  scale_color_gradientn(colors=c("grey","turquoise","blue","navyblue","black")) +
#  
#  # custom axes
#  scale_x_continuous(label = axisdf_splice$chrom,
#                     breaks = axisdf_splice$center, expand=c(.01,.01)) +
#  scale_y_continuous(expand = c(.05,.05)) +
#  
#  # custom theme
#  theme_bw() +
#  theme(legend.position = "none",
#        #axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),
#        text = element_text(size=24),
#        axis.line.x = element_blank(),
#        axis.ticks = element_blank(),
#        axis.text.x=element_blank(),
#        plot.title = element_text(size=18),
#        panel.border = element_blank(),
#        panel.grid.major.x = element_blank(),
#        panel.grid.minor.x = element_blank()) +
#labs(x="")#~paste(Delta,"PSI"))
#splice_manhattan_plot_alt

#ggsave("figures/rMATS_manhattan.png", plot = splice_manhattan_plot,
#       device = "png", width = 8, height = 6, units = "in", dpi = 300)
       

#FIG_2 <- Fst_manhattan_plot / expr_manhattan_plot / splice_manhattan_plot
FIG_2 <- (Fst_manhattan_plot / expr_manhattan_plot / splice_manhattan_plot) | (plot_LFC_vs_Fst / plot_dPSI_vs_Fst)
FIG_2 <- FIG_2 + plot_layout(widths = c(3, 1))
FIG_2

#ggsave("figures/FIG_2_raw.png", plot = FIG_2,
#       device = "png", width = 12, height = 7.2, units = "in", dpi = 300)

ggsave("figures/FIG_2_raw.png", plot = FIG_2,
       device = "png", width = 15, height = 7.2, units = "in", dpi = 300)

#### match with Ath gene names and compare to results from other splicing analyses ####
rmats_ds_genes <- subset(all_AS_events, FDR < .05) %>%
  dplyr::select(GeneID) %>% unique() %>% #just get the unique genes (some genes have multiple events)
  rename(Ha412_gene="GeneID") %>%
  #mutate(Ha412_gene=gsub("gene","mRNA",Ha412_gene)) %>%
  left_join(Ha412_Ath_mappings)

write.table(rmats_ds_genes, file="analysis/rMATS/results_2022-07-14/significant_ds_genes.txt",
            quote = F, row.names = F, sep = "\t")

write.table(rmats_ds_genes$Ha412_gene,
            file = "analysis/GO_analysis/study_DS_rMATS_genes.txt", quote = F,
            row.names = F, sep = "\t", col.names = F)

rmats_dexseq_overlap_genes <- inner_join(rmats_ds_genes, 
                                         deu_genes_noLFCthreshold)
rmats_dexseq_union_genes <- full_join(rmats_ds_genes, deu_genes_noLFCthreshold)
rmats_dexseq_deseq2_overlap_genes <- inner_join(rmats_dexseq_overlap_genes,
                                                de_genes_noLFCthreshold )

parents_diff_ds_genes <- read.table("analysis/GO_analysis/study_DS_parents_diff_genes.txt",
                                    col.names = "Ha412_gene")
all_ds_overlap_genes <- inner_join(rmats_ds_genes, deu_genes_noLFCthreshold) %>%
  inner_join(parents_diff_ds_genes)


write.table(rmats_dexseq_overlap_genes,
            file="analysis/rmats_dexseq_overlap_genes.txt", quote = F,
            row.names = F, sep = "\t")
write.table(rmats_dexseq_overlap_genes$Ha412_gene,
            file="analysis/GO_analysis/study_DS_rMATS_DEU_overlap_genes.txt",
            quote = F, row.names = F, col.names = F)

write.table(rmats_dexseq_union_genes$Ha412_gene,
            file="analysis/GO_analysis/study_DS_rMATS_DEU_union_genes.txt",
            quote = F, row.names = F, col.names = F)

write.table(rmats_dexseq_deseq2_overlap_genes$Ha412_gene, file="analysis/GO_analysis/study_overlap_rMATS_DEU_DE_genes.txt", quote = F, row.names = F, col.names = F)

#### compare fraction of DE–DS overlap genes that are IR events, ####
# compared to the number of IR events in non-overlapping DS genes. 
tmp <- inner_join(dplyr::select(sig_AS_events, ID, Ha412_gene),
                  dplyr::select(de_results_p05_Shrink, Ha412_gene))

tmp <- separate(tmp, ID, into = c("EventType", "ID"), sep="_")

tmp <- tmp %>%
  group_by(EventType) %>%
  summarise(count=n()) %>%
  ggplot(aes(x=EventType,y=count)) +
  geom_bar(stat="Identity")

sum(tmp$count)

tmp2 <- rmatsout_summary %>% subset(Category=="SignificantEventsJCEC")
sum(tmp2$Count)
tmp2

#### plot PSI values for individual events ####
dune_fq <- c("Kane-602-10_S8_R1_001.trimmed.fq.gz",
  "Kane-602-11_S9_R1_001.trimmed.fq.gz",
  "Kane-602-12_S10_R1_001.trimmed.fq.gz",
  "Kane-602-14_S11_R1_001.trimmed.fq.gz",
  "Kane-602-15_S12_R1_001.trimmed.fq.gz",
  "Kane-602-1_S1_R1_001.trimmed.fq.gz",
  "Kane-602-2_S2_R1_001.trimmed.fq.gz",
  "Kane-602-3_S3_R1_001.trimmed.fq.gz",
  "Kane-602-4_S4_R1_001.trimmed.fq.gz",
  "Kane-602-7_S5_R1_001.trimmed.fq.gz",
  "Kane-602-8_S6_R1_001.trimmed.fq.gz",
  "Kane-602-9_S7_R1_001.trimmed.fq.gz")
tmp <- all_AS_events_PSI %>%
  pivot_longer(4:27, names_to = "sample", values_to = "PSI") %>%
  mutate(ecotype=if_else(sample %in% dune_fq, "dune", "non-dune"))

tmp$PSI <- as.numeric(tmp$PSI)

ggplot(data=tmp %>% filter(ID=="RI_2786"),
       aes(x=ecotype, y=PSI,
           fill=ecotype, shape=ecotype)) +
  geom_point(size=5) +
  theme_bw() +
  theme(text = element_text(size=24),#,
        legend.position = "none",
        plot.title = element_text(size=18),
        axis.text = element_text(size=18)) +
  scale_fill_manual(values = c("gold2", "forestgreen")) +
  scale_shape_manual(values=c(24,21))


#### looking at the genes in chr11 peak ####
#chr11_deltaPSI <- all_AS_events_deltaPSI %>% 
#  subset(chr=="chrHa412HOChr11") %>% 
#  subset(FDR < .05) %>% 
#  #subset(abs(IncLevelDifference) > 0.15) %>%
#  dplyr::select(ID, GeneID, IncLevelDifference, FDR)
#
#chr11_deltaPSI
#
## Ha412_Ath_mappings is from merge_gene_names_and_GO_terms.R
#chr11_deltaPSI <- left_join(chr11_deltaPSI, Ha412_Ath_mappings) %>% 
#  dplyr::select(Ath_gene)
#write.table(chr11_deltaPSI, file="analysis/rMATS/chr11_deltaPSI.txt", quote=F, row.names = F)
#