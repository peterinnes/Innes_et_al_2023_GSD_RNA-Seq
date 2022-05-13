# analyzing rMATS output #

#### check the results summary ####
# How many significant DS events of the 5 different types of splicing?
rmatsout_summary <- read.table('analysis/rMATS/summary.txt', header = T)
rmatsout_summary <- pivot_longer(rmatsout_summary,2:9, names_to = "Category", values_to = "Count") %>%
  subset(Category %in% c("SignificantEventsJC", "SignificantEventsJCEC"))

#### merge all the results files together ####
# add splice-type prefix to all the event IDs so we can tell different event types apart after we merge them all together
SE_events <- read.table('analysis/rMATS/SE.MATS.JCEC.txt', header = T) %>%
  mutate(ID = gsub("^", "SE_", ID))
A5SS_events <- read.table('analysis/rMATS/A5SS.MATS.JCEC.txt', header = T) %>%
  mutate(ID = gsub("^", "A5SS_", ID))
A3SS_events <- read.table('analysis/rMATS/A3SS.MATS.JCEC.txt', header = T) %>% 
  mutate(ID = gsub("^", "A3SS_", ID))
MXE_events <- read.table('analysis/rMATS/MXE.MATS.JCEC.txt', header = T) %>%
  mutate(ID = gsub("^", "MXE_", ID))
RI_events <- read.table('analysis/rMATS/RI.MATS.JCEC.txt', header = T) %>% 
  mutate(ID = gsub("^", "RI_", ID))

# count number of significant intron retention splicing events, with a deltaPSI threshold of 15%
RI_events |>
  subset(FDR < .05) |>
  subset(abs(IncLevelDifference) > .15) |>
  dim() 

# combine all the results files into one file, for use in PCA of the Inclusion Levels, in the script pca.R
all_AS_events <- full_join(SE_events, RI_events) %>%
  full_join(MXE_events) %>%
  full_join(A5SS_events) %>% 
  full_join(A3SS_events)

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

# all events data frame should have same number of rows as the sum of all the different event counts in the summary file. And they do. 95,933
dim(all_AS_events)
sum(rmatsout_summary$TotalEventsJCEC)

#### plot deltaPSI across the genome ####
all_AS_events_deltaPSI <- all_AS_events %>%
  mutate(position = coalesce(exonStart_0base, X1stExonStart_0base, longExonStart_0base, riExonStart_0base)) %>%
  dplyr::select(ID, GeneID, chr, position, IncLevelDifference, FDR)

# get the cumulative position of each gene and then the max gene end position on each chrom
cumulative_splice <- all_AS_events_deltaPSI %>% 
  filter(!grepl("Chr00", chr)) %>%
  group_by(chr) %>% 
  summarise(max_bp=max(position)) %>% 
  mutate(bp_add = lag(cumsum(as.numeric(max_bp)), default = 0)) %>% 
  dplyr::select(chr, bp_add)

splice_manhattan_df <- inner_join(all_AS_events_deltaPSI, cumulative_splice) %>%
  rowwise() %>%
  mutate(cumulative_position = position + bp_add, sig=if_else(FDR<.05, "yes", "no"))

axis_set_splice <- splice_manhattan_df %>%
  group_by(chr) %>%
  summarize(max=max(cumulative_position), min=min(cumulative_position)) %>%
  mutate(center=(min+max)/2)

axis_set_splice$chr <- c("01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17")

# split apart significant DE genes and nonsignifcant DE genes so we can downsample the nonsig data, for prettier plot
sig_splice_data <- splice_manhattan_df %>% 
  subset(FDR < 0.05)
notsig_splice_data <- splice_manhattan_df %>% 
  subset(FDR >= 0.05) %>%
  group_by(chr) %>% 
  sample_frac(0.1)

splice_manhattan_df_reduced <- bind_rows(sig_splice_data, notsig_splice_data) 

splice_manhattan_plot <- ggplot(splice_manhattan_df_reduced, aes(x=cumulative_position, y=IncLevelDifference, color=as.factor(chr), alpha=as.factor(sig)=="yes")) +
  geom_point(size=3) +
  scale_x_continuous(label = axis_set_splice$chr, breaks = axis_set_splice$center, expand=c(0.0175,0.0175)) +
  scale_color_manual(values = rep(c("cornflowerblue", "grey50"), unique(length(axis_set_splice$chr)))) +
  scale_alpha_manual(values = c(.1,1)) +
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(size = 36)) +
  geom_hline(yintercept = -.15, linetype = "dashed") +
  geom_hline(yintercept = .15, linetype = "dashed") +
  labs(x="")

# looking at the genes in chr11 peak
chr11_deltaPSI <- all_AS_events_deltaPSI %>% 
  subset(chr=="chrHa412HOChr17") %>% 
  subset(FDR < .05) %>% 
  subset(abs(IncLevelDifference) > 0.15) %>%
  dplyr::select(ID, GeneID, IncLevelDifference, FDR)

head(chr11_deltaPSI)

# Ha412_Ath_mappings is from merge_gene_names_and_GO_terms.R
Ha412_Ath_mappings <- Ha412_Ath_mappings %>%
  mutate(GeneID=gsub("mRNA", "gene", Ha412_gene)) %>%
  dplyr::select(GeneID, Ath_gene)

chr11_deltaPSI <- left_join(chr11_deltaPSI, Ha412_Ath_mappings) %>% 
  dplyr::select(Ath_gene)
write.table(chr11_deltaPSI, file="analysis/rMATS/chr11_deltaPSI.txt", quote=F, row.names = F)
