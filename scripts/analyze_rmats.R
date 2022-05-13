#### analyzing rMATS output ####

# How many significant DS events of the 5 different types of splicing?
rmatsout_summary <- read.table('analysis/rMATS/summary.txt', header = T)
rmatsout_summary <- pivot_longer(rmatsout_summary,2:9, names_to = "Category", values_to = "Count") %>%
  subset(Category %in% c("SignificantEventsJC", "SignificantEventsJCEC"))


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

RI_events |>
  subset(FDR < .05) |>
  subset(abs(IncLevelDifference) > .15) |>
  dim() #matches the counts in summary.txt

all_AS_events <- full_join(SE_events, RI_events) %>%
  full_join(MXE_events) %>%
  full_join(A5SS_events) %>% 
  full_join(A3SS_events) %>%
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
