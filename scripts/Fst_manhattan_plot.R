library(tidyr)

#### make a manhattan plot of genome-wide Fst ####
Fst <- read.table("analysis/Fst/Fst_1Kb_windows.weir.txt", header=T) %>%
  na.omit()

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
cumulative_Fst <- Fst %>% 
  filter(!grepl("Chr00", CHROM)) %>% #filter out non-chromosome contigs
  group_by(CHROM) %>% 
  summarise(max_POS = max(BIN_END)) %>% 
  mutate(POS_add = lag(cumsum(as.numeric(max_POS)), default = 0)) %>% 
  dplyr::select(CHROM, POS_add)

Fst <- Fst %>% 
  inner_join(cumulative, by = "CHROM") %>% 
  mutate(POS_cum = ((BIN_START + BIN_END)/2) + POS_add) %>%
  filter(!grepl("Chr00", CHROM))
#Fst$CHROM <- as.factor(Fst$CHROM)

# get mid point of each chromosome
axis_set_Fst <- Fst %>%
  group_by(CHROM) %>%
  summarize(max=max(POS_cum), min=min(POS_cum)) %>%
  filter(!grepl("Chr00", CHROM)) %>%
  mutate(center=(min+max)/2)
axis_set$CHROM <- c("01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17")

Fst_manhattan_plot <- ggplot(Fst, aes(x=POS_cum, y=MEAN_FST, color=as.factor(CHROM))) +
  geom_point(alpha=0.75) +
  theme_classic() +
  scale_x_continuous(label = axis_set_Fst$CHROM, breaks = axis_set_Fst$center) +
  scale_y_continuous(expand=c(0,0)) +
  scale_color_manual(values = rep(c("cornflowerblue", "grey50"), unique(length(axis_set_Fst$CHROM)))) +
  coord_cartesian(clip = 'off') +
  labs(x="Chromosome", y=expression(italic(F)["st"])) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        text = element_text(size=36),
        axis.line.x = element_blank(),
        panel.border = element_blank(),
        plot.margin = unit(c(.5,0,.25,.25), 'cm'))
  
ggsave("figures/Fst_manhattan.jpg", plot=Fst_manhattan_plot, width=28, height = 14, dpi = 600, units = "cm")
