# analyzing rMATS output #
library(tidyr)
library(ggplot2)

#### check the results summary ####
# How many significant DS events of the 5 different types of splicing?
rmatsout_summary <- read.table('analysis/rMATS/results_2022-07-14/summary.txt', header = T)
rmatsout_summary <- pivot_longer(rmatsout_summary,2:9, names_to = "Category", values_to = "Count") %>%
  subset(Category %in% c("TotalEventsJCEC", "SignificantEventsJCEC"))
ggplot(data = rmatsout_summary %>% subset(Category=="SignificantEventsJCEC"),
       aes(x=EventType, y=Count)) +
  geom_bar(stat = 'identity')
  

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
  subset(abs(IncLevelDifference) > .05) |>
  dim() 

# combine all the results files into one file, for use in PCA of the Inclusion Levels, in the script pca.R
all_AS_events <- full_join(SE_events, RI_events) %>%
  full_join(MXE_events) %>%
  full_join(A5SS_events) %>% 
  full_join(A3SS_events)

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
non.dune_fq <- c("Kane-602-16_S13_R1_001.trimmed.fq.gz",
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
                 "Kane-602-30_S24_R1_001.trimmed.fq.gz")

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

dim(all_AS_events_PSI)
#sum(rmatsout_summary$TotalEventsJCEC)

sig_AS_events <- subset(all_AS_events, FDR<.05)

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

#### plot deltaPSI across the genome ####
all_AS_events_deltaPSI <- all_AS_events %>%
  mutate(exon_start = coalesce(exonStart_0base, X1stExonStart_0base, longExonStart_0base, riExonStart_0base)) %>%
  dplyr::select(ID, GeneID, "chrom"=chr, exon_start, IncLevelDifference, FDR) %>%
  mutate(chrom=gsub("chr", "", chrom)) %>%
  arrange(chrom, exon_start)

# get the cumulative length of each chromosome
chrom_lengths <- read.table("Ha412H0_chrom_lengths.txt", col.names = c("chrom", "length")) %>%
  mutate(length=as.numeric(length))
cum_lengths <- chrom_lengths %>%
  mutate(tot=cumsum(length)-length) %>%
  dplyr::select(-length)

# get inversion regions
inversion_regions <- read.table("analysis/inversions/inversion_regions.txt") %>%
  dplyr::select(V2,V3,V4,V9) %>%
  rename(V2="chrom", V3="start", V4="end", V9="name") %>%
  left_join(cum_lengths) %>%
  mutate(startcum=start+tot, endcum=end+tot)

splice_manhattan_df <- inner_join(all_AS_events_deltaPSI, cum_lengths) %>%
  #rowwise() %>%
  mutate(exon_startcum = exon_start + tot, sig=if_else(FDR<.05, "yes", "no"))

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

splice_manhattan_plot <- ggplot(splice_manhattan_df_reduced, aes(x=exon_startcum, y=IncLevelDifference, color=as.factor(chrom), alpha=as.factor(sig)=="yes")) +
  geom_point(size=3) +
  scale_x_continuous(label = axisdf_splice$chrom, breaks = axisdf_splice$center, expand=c(0.0175,0.0175)) +
  #scale_y_continuous(expand=expansion(mult=c(0,0))) +
  scale_color_manual(values = rep(c("cornflowerblue", "grey50"), unique(length(axisdf_splice$chrom)))) +
  scale_alpha_manual(values = c(.1,1)) +
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(size = 36),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  #geom_hline(yintercept = -.05, linetype = "dashed") +
  #geom_hline(yintercept = .05, linetype = "dashed") +
  geom_segment(data=inversion_regions, aes(x=startcum, xend=endcum, y=0, yend=0),
               size=6, color="red", alpha=0.8) +
  labs(x="Chromosome", title="Differential splicing\n(event-based)")
ggsave("figures/rMATS_manhattan.png", plot = splice_manhattan_plot, device = "png", width = 8, height = 6, units = "in", dpi = 300)
       

#### match with Ath gene names and compare to DEXSeq results ####
rmats_ds_genes <- subset(all_AS_events, FDR < .05) %>%
  dplyr::select(GeneID) %>% unique() %>% #just get the unique genes (some genes have multiple events)
  rename("GeneID"="Ha412_gene") %>%
  #mutate(Ha412_gene=gsub("gene","mRNA",Ha412_gene)) %>%
  left_join(Ha412_Ath_mappings)

write.table(rmats_ds_genes, file="analysis/rMATS/results_2022-07-14/significant_ds_genes.txt", quote = F, row.names = F, sep = "\t")

write.table(rmats_ds_genes$Ha412_gene, file = "analysis/GO_analysis/study_DS_rMATS_genes.txt", quote = F, row.names = F, sep = "\t", col.names = F)

rmats_dexseq_overlap_genes <- inner_join(rmats_ds_genes, 
                                         deu_genes_noLFCthreshold)
rmats_dexseq_union_genes <- full_join(rmats_ds_genes, deu_genes_noLFCthreshold)
rmats_dexseq_deseq2_overlap_genes <- inner_join(rmats_dexseq_overlap_genes, de_genes_noLFCthreshold )

#ds_overlap_genes <- inner_join(rmats_dexseq_overlap_genes, set4) # set4 is parents_diffv2 genes


write.table(rmats_dexseq_overlap_genes, file="analysis/rmats_dexseq_overlap_genes.txt", quote = F, row.names = F, sep = "\t")
write.table(rmats_dexseq_overlap_genes$Ha412_gene, file="analysis/GO_analysis/study_DS_rMATS_DEU_overlap_genes.txt", quote = F, row.names = F, col.names = F)

write.table(rmats_dexseq_union_genes$Ha412_gene, file="analysis/GO_analysis/study_DS_rMATS_DEU_union_genes.txt", quote = F, row.names = F, col.names = F)

write.table(rmats_dexseq_deseq2_overlap_genes$Ha412_gene, file="analysis/GO_analysis/study_overlap_rMATS_DEU_DE_genes.txt", quote = F, row.names = F, col.names = F)

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