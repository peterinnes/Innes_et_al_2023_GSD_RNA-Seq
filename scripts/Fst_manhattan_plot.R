library(tidyr)

#### make a manhattan plot of genome-wide Fst ####
Fst <- read.table("analysis/Fst/Fst_1Kb_windows.weir.txt", header=T) %>%
  na.omit() %>%
  rename(CHROM="chrom")
Fst <- read.table("analysis/Fst/Fst_10Kb_windows.weir.txt", header=T) %>%
  na.omit() %>%
  rename(CHROM="chrom")

# convert negative weighted Fst values to 0. Okay to do this?
for (i in 1:length(Fst$WEIGHTED_FST)) {
  if (Fst$WEIGHTED_FST[i] < 0) {
   Fst$WEIGHTED_FST[i] <- 0
  }
  if (Fst$MEAN_FST[i] < 0) {
     Fst$MEAN_FST[i] <- 0
  }
}



# get the cumulative position of each SNP/window and then the max SNP end position on each chrom
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

Fst_manhattan_df <- Fst %>% 
  inner_join(cum_lengths, by = "chrom") %>% 
  mutate(cumstart = BIN_START + tot) %>%
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

Fst_manhattan_plot <- ggplot(Fst_manhattan_df, aes(x=cumstart, y=MEAN_FST, color=as.factor(chrom))) +
  geom_point(alpha=0.75, size=3) +
  theme_classic() +
  scale_x_continuous(label = axisdf_Fst$chrom, breaks = axisdf_Fst$center) +
  scale_y_continuous(expand=c(0,0)) +
  scale_color_manual(values = rep(c("cornflowerblue", "grey50"), unique(length(axisdf_Fst$chrom)))) +
  coord_cartesian(clip = 'off') +
  labs(x="", y=expression(italic(F)["st"]), title="Genetic divergence") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),
        text = element_text(size=36),
        axis.line.x = element_blank(),
        panel.border = element_blank()) + #,
        #plot.margin = unit(c(.5,0,.25,.25), 'in')) +
  # add inversion regions
  geom_segment(data=inversion_regions, aes(x=startcum, xend=endcum, y=0, yend=0),
               size=6, color="red", alpha=0.8)
  
ggsave("figures/Fst_manhattan.png", device="png",plot=Fst_manhattan_plot, width=8, height = 6, dpi = 300, units = "in")
