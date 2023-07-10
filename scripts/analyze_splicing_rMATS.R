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

# plot non-significant events and significant events together. for sup mat
plot_rmats_event_counts <- ggplot(data = rmatsout_summary %>% subset(Category=="TotalEventsJCEC"),
       aes(x=EventType, y=Count)) +
  geom_col(position = 'identity', color="grey50", alpha=0.5) +
  geom_col(data=rmatsout_summary %>% subset(Category=="SignificantEventsJCEC"),
           aes(x=EventType, y=Count),fill="black", position = 'identity') +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  labs(x="Event Type")

ggsave("figures/plot_rMATS_event_counts_raw.pdf", plot=plot_rmats_event_counts,
       device = "pdf", height = 87.5, width = 87.5, dpi = 300, units = "mm")

# significant events bar plot
ggplot(data = rmatsout_summary %>% subset(Category=="SignificantEventsJCEC"),
       aes(x=EventType, y=Count)) +
  geom_bar(stat = 'identity') +
  theme_classic() +
  theme(text = element_text(size=24),
        axis.text = element_text(size=18))

#### merge all the results files together ####
# add splice-type prefix to all the event IDs so we can tell different event types apart after we merge them all together
SE_events <- read.table('analysis/rMATS/results_2022-07-14/SE.MATS.JCEC.txt',
                        header = T) %>% mutate(ID = gsub("^", "SE_", ID))
A5SS_events <- read.table('analysis/rMATS/results_2022-07-14/A5SS.MATS.JCEC.txt',
                          header = T) %>% mutate(ID = gsub("^", "A5SS_", ID))
A3SS_events <- read.table('analysis/rMATS/results_2022-07-14/A3SS.MATS.JCEC.txt',
                          header = T) %>% mutate(ID = gsub("^", "A3SS_", ID))
MXE_events <- read.table('analysis/rMATS/results_2022-07-14/MXE.MATS.JCEC.txt',
                         header = T) %>% mutate(ID = gsub("^", "MXE_", ID))
RI_events <- read.table('analysis/rMATS/results_2022-07-14/RI.MATS.JCEC.txt',
                        header = T) %>% mutate(ID = gsub("^", "RI_", ID))

# count number of significant intron retention splicing events, Should match the summary
RI_events |>
  subset(FDR < .05) |>
  #subset(abs(IncLevelDifference) > .05) |>
  dim() 

## attempting to look at novel Splice Sites
#RI_novelSS <- read.table('analysis/rMATS/results_2022-07-14/fromGTF.novelSpliceSite.RI.txt',
#                         header = T) %>% mutate(ID = gsub("^", "RI_", ID))
#SE_novelSS <- read.table('analysis/rMATS/results_2022-07-14/fromGTF.novelSpliceSite.SE.txt',
#                         header = T) %>% mutate(ID = gsub("^", "SE_", ID))
#A5SS_novelSS <- read.table('analysis/rMATS/results_2022-07-14/fromGTF.novelSpliceSite.A5SS.txt',
#                           header = T) %>% mutate(ID = gsub("^", "A5SS_", ID))
#A3SS_novelSS <- read.table('analysis/rMATS/results_2022-07-14/fromGTF.novelSpliceSite.A3SS.txt',
#                           header = T) %>% mutate(ID = gsub("^", "A3SS_", ID))
#MXE_novelSS <- read.table('analysis/rMATS/results_2022-07-14/fromGTF.novelSpliceSite.MXE.txt',
#                           header = T) %>% mutate(ID = gsub("^", "MXE_", ID))
#
## fraction of events that are based on novel splice sites
#inner_join(SE_events, SE_novelSS) %>% nrow() / nrow(SE_events)
#inner_join(A5SS_events, A5SS_novelSS) %>% nrow() / nrow(A5SS_events)
#inner_join(A3SS_events, A3SS_novelSS) %>% nrow() / nrow(A3SS_events)
#inner_join(RI_events, RI_novelSS) %>% nrow() / nrow(RI_events)
#inner_join(MXE_events, MXE_novelSS) %>% nrow() / nrow(MXE_events)
#
## get postiions of novel splice site regions, following SnpEff +8/-3 range for
## definiing a splice site ¨region¨ i.e. a small window encompassing a splice site. 
## if plus strand, "region 1" is the acceptor; "region 2" is donor
## if minus strand, "region 2" is the acceptor; "region 1" is donor
#SE_novel_splice_regions <- inner_join(SE_events, SE_novelSS) %>%
#  dplyr::select(ID, chr, GeneID, strand, exonStart_0base, exonEnd) %>%
#  mutate(splice_region_1_start = exonStart_0base - 8,
#         splice_region_1_end = exonStart_0base + 3,
#         splice_region_2_start = exonEnd - 3,
#         splice_region_2_end = exonEnd + 8)
#
#A5SS_novel_splice_regions <- inner_join(A5SS_events, A5SS_novelSS) %>% 
#  dplyr::select(ID, chr, GeneID, strand, longExonStart_0base, longExonEnd) %>%
#  mutate(splice_region_1_start = longExonStart_0base - 8,
#         splice_region_1_end = longExonStart_0base + 3,
#         splice_region_2_start = longExonEnd - 3,
#         splice_region_2_end = longExonEnd + 8)
#
#A3SS_novel_splice_regions <- inner_join(A3SS_events, A3SS_novelSS) %>% 
#  dplyr::select(ID, chr, GeneID, strand, longExonStart_0base, longExonEnd) %>%
#  mutate(splice_region_1_start = longExonStart_0base - 8,
#         splice_region_1_end = longExonStart_0base + 3,
#         splice_region_2_start = longExonEnd - 3,
#         splice_region_2_end = longExonEnd + 8)
#
#MXE_novel_splice_regions <- inner_join(MXE_events, MXE_novelSS) %>%
#  dplyr::select(ID, chr, GeneID, strand,
#                X1stExonStart_0base, X1stExonEnd,
#                X2ndExonStart_0base, X2ndExonEnd) %>%
#  mutate(splice_region_1_start = X1stExonStart_0base - 8,
#         splice_region_1_end = X1stExonStart_0base + 3,
#         splice_region_2_start = X1stExonEnd - 3,
#         splice_region_2_end = X1stExonEnd + 8,
#         splice_region_3_start = X2ndExonStart_0base - 8,
#         splice_region_3_end = X2ndExonStart_0base + 3,
#         splice_region_4_start = X2ndExonEnd - 3,
#         splice_region_4_end = X2ndExonEnd + 8)
#
## for intron retention, mutations within the retained intron might cause RI
## set the window to be the intron that is retained.
#RI_novel_splice_regions <- inner_join(RI_events, RI_novelSS) %>%
#  dplyr::select(ID, chr, GeneID, strand, upstreamEE, downstreamES) %>%
#  mutate(splice_region_1_start=upstreamEE-3,
#         splice_region_1_end=downstreamES+3)
#
## combine all novel splice regions
#novel_splice_regions <- full_join(SE_novel_splice_regions,
#                                  RI_novel_splice_regions) %>%
#  full_join(MXE_novel_splice_regions) %>%
#  full_join(A5SS_novel_splice_regions) %>%
#  full_join(A3SS_novel_splice_regions) %>%
#  dplyr::select(chr, GeneID, ID, strand,
#                splice_region_1_start, splice_region_1_end,
#                splice_region_2_start, splice_region_2_end,
#                splice_region_3_start, splice_region_3_end,
#                splice_region_4_start, splice_region_4_end)
#
#region_1 <- novel_splice_regions %>%
#  dplyr::select(chr, GeneID, ID, splice_region_1_start, splice_region_1_end) %>%
#  dplyr::rename(start=splice_region_1_start, end=splice_region_1_end) %>%
#  mutate(region="region_1")
#
#region_2 <- novel_splice_regions %>%
#  dplyr::select(chr, GeneID, ID, splice_region_2_start, splice_region_2_end) %>%
#  dplyr::rename(start=splice_region_2_start, end=splice_region_2_end) %>%
#  mutate(region="region_2") %>%
#  na.omit()
#
#region_3 <- novel_splice_regions %>%
#  dplyr::select(chr, GeneID, ID, splice_region_3_start, splice_region_3_end) %>%
#  dplyr::rename(start=splice_region_3_start, end=splice_region_3_end) %>%
#  mutate(region="region_3") %>%
#  na.omit()
#
#region_4 <- novel_splice_regions %>%
#  dplyr::select(chr, GeneID, ID, splice_region_4_start, splice_region_4_end) %>%
#  dplyr::rename(start=splice_region_4_start, end=splice_region_4_end) %>%
#  mutate(region="region_4") %>%
#  na.omit()
#
#novel_splice_regions <- full_join(region_1, region_2) %>%
#  full_join(region_3) %>%
#  full_join(region_4) %>% dplyr::select(chr, start, end, GeneID, EventID=ID, region)
#
#write.table(novel_splice_regions,
#            file = "analysis/snpEff/rMATS_novel_splice_regions.txt",
#            sep = "\t", quote = F, row.names = F)

# combine all the events results files into one file,
all_AS_events <- full_join(SE_events, RI_events) %>%
  full_join(MXE_events) %>%
  full_join(A5SS_events) %>% 
  full_join(A3SS_events) %>%
  arrange(FDR) %>%
  left_join(Ha412_Ath_mappings %>% dplyr::rename(GeneID=Ha412_gene)) %>%
  relocate(ID,GeneID, Ath_gene)

# split apart valules for individual plants in the IncLevel (PSI) columns 
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

# filter results for missingness (remove events with > 40% missingness)
# for use in PCA of the Inclusion Levels, in the script pca.R
all_AS_events_PSI[,4:27] <- all_AS_events_PSI[,4:27] %>%
  lapply(as.numeric) %>% as.data.frame()
tmp <- all_AS_events_PSI[,4:27]
rownames(tmp) <- all_AS_events_PSI$ID
rmats_splice_df_filtered <- tmp %>%
  mutate(missing_perc = rowMeans(is.na(tmp))) %>%
  filter(missing_perc<.4) %>%
  dplyr::select(!missing_perc)

# set remaining missing values to the average PSI of each event (each row)
for( i in 1:nrow(rmats_splice_df_filtered) ){
  for( j in 1:ncol(rmats_splice_df_filtered) ){
    if( is.na(rmats_splice_df_filtered[i,j]) ){
      rmats_splice_df_filtered[i,j] <- rowMeans(rmats_splice_df_filtered[i,],
                                                na.rm = T)
    }
  }
}

#get_log2FC <- function(A,B){
#  log2FC <- log2(A) - log2(B)
#  return(log2FC)
#}
## calculate log2FoldChange for dune non-dune 
#all_AS_events_PSI$log2FC <- get_log2FC(rowMeans(all_AS_events_PSI[4:15], na.rm = T),
#                                       rowMeans(all_AS_events_PSI[16:27], na.rm = T))

all_AS_events_deltaPSI <- all_AS_events %>%
  mutate(event_start = coalesce(exonStart_0base, X1stExonStart_0base,
                               longExonStart_0base, upstreamEE),
         event_end = coalesce(exonEnd, X2ndExonEnd,
                             longExonEnd, downstreamES)) %>%
  dplyr::select(ID, GeneID,Ath_gene, "chrom"=chr, event_start, event_end, IncLevelDifference, FDR) %>%
  mutate(chrom=gsub("chr", "", chrom)) %>%
  filter(ID %in% rownames(rmats_splice_df_filtered)) #filter out events with > 40% missingness

# save as Rdata for use in PCA and manhattan plot
save(rmats_splice_df_filtered, all_AS_events_PSI, all_AS_events_deltaPSI, 
     file = "data2/Rdata/rmats_results_dfs.Rdata")

# save to file while filtering for < 40% missingness
write.table(all_AS_events %>% filter(ID %in% rownames(rmats_splice_df_filtered)),
            "analysis/rMATS/results_2022-07-14/all_AS_events.txt", sep = '\t',
            quote = F, row.names = F)
write.table(all_AS_events_PSI %>% filter(ID %in% rownames(rmats_splice_df_filtered)),
            "analysis/rMATS/results_2022-07-14/all_AS_events_PSI.txt",
            row.names = F, quote = F, sep = '\t')
# save to file
write.table(all_AS_events_deltaPSI,
            "analysis/rMATS/results_2022-07-14/all_AS_events_deltaPSI.tsv",
            sep = '\t', quote = F, row.names = F)

#for( i in 1:nrow(all_AS_events_deltaPSI)){
#  if( all_AS_events_deltaPSI$FDR[i] == 0 ){
#    all_AS_events_deltaPSI$FDR[i] <- 1.134807e-14 #next smallest FDR
#  }
#}


#### does the dune ecotype retain more introns? ####
RI_dPSI <- filter(all_AS_events_deltaPSI, grepl("RI",ID) & FDR < .05)
dim(RI_dPSI %>% filter(IncLevelDifference > 0))
dim(RI_dPSI %>% filter(IncLevelDifference < 0))

SE_dPSI <- filter(all_AS_events_deltaPSI, grepl("SE",ID) & FDR < .05)
dim(SE_dPSI %>% filter(IncLevelDifference > 0))
dim(SE_dPSI %>% filter(IncLevelDifference < 0))

#### splicing manhattan plot ####
# get the cumulative length of each chromosome
chrom_lengths <- read.table("data/ref_genome_Ha412HO/chrom_sizes_Ha412HOv2.0.txt",
                            col.names = c("chrom", "length")) %>%
  mutate(length=as.numeric(length))

cum_lengths <- chrom_lengths %>%
  mutate(cumstart=cumsum(length)-length,
         cumend=cumsum(length))

# get inversion regions
inversion_regions <- read.table("analysis/inversions/inversion_regions.txt") %>%
  dplyr::select(V2,V3,V4,V9)
names(inversion_regions) <- c("chrom", "start", "end", "name")
inversion_regions <- inversion_regions %>%
  filter(name %in% c("pet05.01", "pet09.01", "pet11.01", "pet17.01")) %>%
  left_join(cum_lengths[c(1,3)]) %>%
  mutate(inv_cumstart=start+cumstart, inv_cumend=end+cumstart)

splice_manhattan_df <- inner_join(all_AS_events_deltaPSI, cum_lengths) %>%
  #rowwise() %>%
  mutate(event_startcum = event_start + cumstart, sig=if_else(FDR<.05, "yes", "no"))

axisdf_splice <- splice_manhattan_df %>%
  group_by(chrom) %>%
  summarize(center=(max(event_startcum) + min(event_startcum)) / 2) %>%
  mutate(chrom=c("1","2","3","4","5",
                 "6","7","8","9","10",
                 "11","12","13","14","15",
                 "16","17"))

# split apart significant DE genes and nonsignifcant DE genes so we can downsample the nonsig data, for better vis
sig_splice_data <- splice_manhattan_df %>% 
  subset(FDR < 0.05)
#sig_splice_data <- splice_manhattan_df %>% 
#  subset(FDR < 0.05) %>%
#  group_by(GeneID) %>%
#  slice_max(abs(IncLevelDifference))
notsig_splice_data <- splice_manhattan_df %>% 
  subset(FDR >= 0.05) %>%
  group_by(chrom) %>% 
  sample_frac(0.5)

splice_manhattan_df_reduced <- bind_rows(sig_splice_data, notsig_splice_data) 

splice_manhattan_plot <- ggplot(splice_manhattan_df_reduced,
                                aes(x=event_startcum, y=IncLevelDifference),
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
  scale_size_manual(name="Significance (FDR < .05)", values = c(2,.5), labels=c("yes", "no")) +
  #geom_point(alpha=0.75, aes(size=-log(FDR))) +
  #scale_color_gradientn(colors=c("grey","turquoise","blue","navyblue","black")) +
  
  # custom axes
  scale_x_continuous(label = axisdf_splice$chrom,
                     breaks = axisdf_splice$center, expand=c(0.01,0.01)) +
  scale_y_continuous(expand=c(0.01,0.01)) +
  
  # custom theme
  theme_bw(base_size = 18) +
  theme(legend.position = "top",
        #axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),
        #text = element_text(size=24),
        axis.line.x = element_blank(),
        axis.ticks = element_blank(),
        #axis.text.x=element_blank(),
        #plot.title = element_text(size=18),
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.text = element_text(size=14),
        legend.title = element_text(size=14),
        #axis.text = element_text(size=18),
        legend.box.margin = margin(t = -30)) +
  
  labs(x="Chromosome",y=~paste(Delta,"PSI"))
splice_manhattan_plot

# without the scatterplots (panels D and E)
FIG_2 <- Fst_manhattan_plot / expr_manhattan_plot / splice_manhattan_plot

# with the scatterplots (panels D and E)
FIG_2 <- (Fst_manhattan_plot / expr_manhattan_plot / splice_manhattan_plot) / (plot_spacer() |plot_LFC_vs_Fst | plot_dPSI_vs_Fst | plot_spacer()) + plot_layout(heights = c(1,1,1,1.25))

FIG_2
FIG_2 <- FIG_2 + plot_layout(widths = c(3, 1))


ggsave("figures/FIG_2_raw.png", plot = FIG_2,
       device = "png", width = 12, height = 10, units = "in", dpi = 300)

# need to revise the plot font size and maybe patchwork widths in order to 
# have the plot look right at the required dimensions (max 175mm).
ggsave("figures/FIG_2_raw.pdf", plot = FIG_2,
       device = "pdf", width = 175, height = 145.833, units = "mm", dpi = 300)

ggsave("figures/plot_splice_manhattan.pdf", plot = splice_manhattan_plot,
       device = "pdf", width = 175, units = "mm", dpi = 300)


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
PSI_df <- all_AS_events_PSI %>%
  pivot_longer(4:27, names_to = "sample", values_to = "PSI") %>%
  mutate(ecotype=if_else(sample %in% dune_fq, "Dune", "Non-dune"))

PSI_df$PSI <- as.numeric(PSI_df$PSI)

# GLH17
plot_GLH17_psi <- ggplot(data=PSI_df %>% filter(ID=="SE_8901"),
       aes(x=ecotype, y=PSI,
           fill=ecotype, shape=ecotype)) +
  geom_point(size=5, alpha=.75, position = position_jitter(width = .1)) +
  geom_boxplot(alpha=0.25, outlier.shape = NA) +
  theme_bw(base_size = 18) +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("gold2", "forestgreen")) +
  scale_color_manual(values = c("gold2", "forestgreen")) +
  scale_shape_manual(values=c(24,21)) +
  labs(x="Ecotype")

# Ha412HOChr14g0679271. CHLD. TRINITY_DN20639_c0_g1

plot_ALB1_psi <- ggplot(data=PSI_df %>% filter(ID=="RI_21314"),
                         aes(x=ecotype, y=PSI,
                             fill=ecotype, shape=ecotype)) +
  geom_point(size=5, alpha=.75, position = position_jitter(width = .1)) +
  geom_boxplot(alpha=0.25, outlier.shape = NA) +
  theme_bw(base_size = 18) +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("gold2", "forestgreen")) +
  scale_color_manual(values = c("gold2", "forestgreen")) +
  scale_shape_manual(values=c(24,21)) +
  labs(x="Ecotype")

# Ha412HOChr09g0402961. A5SS_13302
plot_ABCK14_psi <- ggplot(data=PSI_df %>% filter(ID=="A5SS_13302"),
                        aes(x=ecotype, y=PSI,
                            fill=ecotype, shape=ecotype)) +
  geom_point(size=5, alpha=.75, position = position_jitter(width = .1)) +
  geom_boxplot(alpha=0.25, outlier.shape = NA) +
  theme_bw(base_size = 18) +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("gold2", "forestgreen")) +
  scale_color_manual(values = c("gold2", "forestgreen")) +
  scale_shape_manual(values=c(24,21)) +
  labs(x="Ecotype")

# rhodanese gene Ha412HOChr10g0445261 A5SS_14127.
plot_Chr10g0445261_psi <- ggplot(data=PSI_df %>% filter(ID=="A5SS_14127"),
                         aes(x=ecotype, y=PSI,
                             fill=ecotype, shape=ecotype)) +
  geom_point(size=5, alpha=.75, position = position_jitter(width = .1)) +
  geom_boxplot(alpha=0.25, outlier.shape = NA) +
  theme_bw(base_size = 18) +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("gold2", "forestgreen")) +
  scale_color_manual(values = c("gold2", "forestgreen")) +
  scale_shape_manual(values=c(24,21)) +
  labs(x="Ecotype")

#Ha412HOChr11g0508211
plot_COR15A_psi <- ggplot(data=PSI_df %>% filter(ID=="SE_17019"),
                          aes(x=ecotype, y=PSI,
                              fill=ecotype, shape=ecotype)) +
  geom_point(size=5, alpha=.75, position = position_jitter(width = .1)) +
  geom_boxplot(alpha=0.25, outlier.shape = NA) +
  theme_bw(base_size = 18) +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("gold2", "forestgreen")) +
  scale_color_manual(values = c("gold2", "forestgreen")) +
  scale_shape_manual(values=c(24,21)) +
  labs(x="Ecotype")

plot_CESA6_psi <- ggplot(data=PSI_df %>% filter(ID=="RI_15992"),
                         aes(x=ecotype, y=PSI,
                             fill=ecotype, shape=ecotype)) +
  geom_point(size=5, alpha=.75, position = position_jitter(width = .1)) +
  geom_violin(alpha=0.25) +
  theme_bw(base_size = 18) +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("gold2", "forestgreen")) +
  scale_color_manual(values = c("gold2", "forestgreen")) +
  scale_shape_manual(values=c(24,21)) +
  labs(x="Ecotype")

plot_EMB2004_psi <- ggplot(data=PSI_df %>% filter(ID=="RI_583"),
                         aes(x=ecotype, y=PSI,
                             fill=ecotype, shape=ecotype)) +
  geom_point(size=5, alpha=.75, position = position_jitter(width = .1)) +
  geom_violin(alpha=0.25) +
  theme_bw(base_size = 18) +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("gold2", "forestgreen")) +
  scale_color_manual(values = c("gold2", "forestgreen")) +
  scale_shape_manual(values=c(24,21)) +
  labs(x="Ecotype")

# Ha412HOChr15g0742211; transcriptome match is TRINITY_DN11324_c1_g1
# RI_22481 is the RI event with highest dPSI for this gene. 
plot_SCPL13_psi <- ggplot(data=PSI_df %>% filter(ID=="RI_22481"),
                         aes(x=ecotype, y=PSI,
                             fill=ecotype, shape=ecotype)) +
  geom_point(size=5, alpha=.75, position = position_jitter(width = .1)) +
  geom_violin(alpha=0.25) +
  theme_bw(base_size = 18) +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("gold2", "forestgreen")) +
  scale_color_manual(values = c("gold2", "forestgreen")) +
  scale_shape_manual(values=c(24,21)) +
  labs(x="Ecotype")
  
#### compare fraction of DE–DS overlap genes that are IR events, ####
# compared to the number of IR events in non-overlapping DS genes.

sig_AS_events <- subset(all_AS_events, FDR<.05) %>%
  rename("Ha412_gene"=GeneID)
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

#### looking at the genes in chr11 inversion ####
chr11_deltaPSI <- all_AS_events_deltaPSI %>% 
  subset(chr=="chrHa412HOChr11") %>% 
  subset(FDR < .05) %>% 
  subset(abs(IncLevelDifference) > 0.15) %>%
  dplyr::select(ID, GeneID, IncLevelDifference, FDR)

#chr11_deltaPSI
#
## Ha412_Ath_mappings is from merge_gene_names_and_GO_terms.R
#chr11_deltaPSI <- left_join(chr11_deltaPSI, Ha412_Ath_mappings) %>% 
#  dplyr::select(Ath_gene)
#write.table(chr11_deltaPSI, file="analysis/rMATS/chr11_deltaPSI.txt", quote=F, row.names = F)
#