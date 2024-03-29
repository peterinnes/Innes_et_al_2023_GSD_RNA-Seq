library(ggplot2)
library(patchwork)
library(dplyr)

#### read-in phenotype data ####
a <- read.csv("~/Desktop/SUNFLOWER/gsd_RNA-seq/GSDgreenhouseExpt2018_Data_sheet1.csv",
              header = T) %>%
  dplyr::select(HABITAT, POPULATION, plant_number, germinated, planted,
                GERM_DATE, ACHENE_LEN, HYPO_LEN, ROOT_LEN, HEIGHT, PLANT_MASS,
                LEAF_num, LEAF_MASS_WET, LEAF_MASS_DRY, RNA_EXT, COMMENTS)
a$LEAF_num[160] <- 8 #fixing data entry error (leaf number entered as 18, should be 8)

b <- read.csv("~/Desktop/SUNFLOWER/gsd_RNA-seq/GSDgreenhouseExpt2018_Data_sheet2.csv",
              header = T) %>%
  dplyr::select(HABITAT, POPULATION, plant_number, code_number, RNA_SEQ)

phenotype_data <- full_join(a,b) %>%
  mutate(HABITAT=gsub("dune","Dune",HABITAT),
         HABITAT=gsub("nonDune","Non-dune",HABITAT)) %>%
  filter(planted==1)

save(phenotype_data, file = "data2/Rdata/seedling_phenotype_data.Rdata")
#### all plants plots ####
height_plot <- ggplot(phenotype_data,
                      aes(x=HABITAT, y=HEIGHT,
                          shape=HABITAT, fill=HABITAT)) +

  geom_point(size=1, alpha=.75, stroke=.25,
             position = position_jitter(width = .2)) +
  geom_violin(alpha=0.25, lwd=.25) +
  theme_bw(base_size = 8) +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("gold2", "forestgreen")) +
  scale_color_manual(values = c("gold2", "forestgreen")) +
  scale_shape_manual(values=c(24,21)) +
  labs(x="", y="Height (cm)")

leaf_num_plot <- ggplot(phenotype_data,
                        aes(x=HABITAT, y=LEAF_num,
                            shape=HABITAT, fill=HABITAT)) +
  geom_point(size=1, alpha=.75, stroke=.25,
             position = position_jitter(width = .2, height = 0)) +
  geom_violin(alpha=0.25, lwd=.25) +
  theme_bw(base_size = 8) +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("gold2", "forestgreen")) +
  scale_color_manual(values = c("gold2", "forestgreen")) +
  scale_shape_manual(values=c(24,21)) +
  labs(x="", y="Num. leaves")

leaf_mass_dry_plot <- ggplot(phenotype_data,
                             aes(x=HABITAT, y=LEAF_MASS_DRY,
                                 shape=HABITAT, fill=HABITAT)) +
  
  geom_point(size=1, alpha=.75, stroke=.25,
             position = position_jitter(width = .2)) +
  geom_violin(alpha=0.25, lwd=.25) +
  theme_bw(base_size = 8) +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("gold2", "forestgreen")) +
  scale_color_manual(values = c("gold2", "forestgreen")) +
  scale_shape_manual(values=c(24,21)) +
  labs(x="Ecotype", y="Leaf dry mass (g)")

plant_mass_plot <- ggplot(phenotype_data,
                        aes(x=HABITAT, y=PLANT_MASS,
                            shape=HABITAT, fill=HABITAT)) +
  geom_point(size=1, alpha=.75, stroke=.25,
             position = position_jitter(width = .2, height = 0)) +
  geom_violin(alpha=0.25, lwd=.25) +
  theme_bw(base_size = 8) +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("gold2", "forestgreen")) +
  scale_color_manual(values = c("gold2", "forestgreen")) +
  scale_shape_manual(values=c(24,21)) +
  labs(x="Ecotype", y="Total dry mass (g)")


#
#leaf_mass_proportion_plot <- ggplot(phenotype_data,
#                          aes(x=HABITAT, y=LEAF_MASS_DRY/PLANT_MASS,
#                              shape=HABITAT, fill=HABITAT)) +
#  geom_point(size=1, alpha=.75, stroke=.25,
#             position = position_jitter(width = .2, height = 0)) +
#  geom_violin(alpha=0.25, lwd=.25) +
#  theme_bw(base_size = 8) +
#  theme(legend.position = "none") +
#  scale_fill_manual(values = c("gold2", "forestgreen")) +
#  scale_color_manual(values = c("gold2", "forestgreen")) +
#  scale_shape_manual(values=c(24,21))

plot_map_and_seedling_traits <- gsd_sampling_map | (height_plot | leaf_num_plot) / (leaf_mass_dry_plot | plant_mass_plot)
ggsave("figures/map_and_seedling_traits.pdf", plot=plot_map_and_seedling_traits,
       device = "pdf", height = 86.625, width = 115.5, dpi = 300, units = "mm")
#seedling_trait_plot <- (height_plot | leaf_num_plot) / (leaf_mass_dry_plot | plant_mass_plot)
#ggsave("figures/seedling_traits_raw.pdf", plot=seedling_trait_plot, device = "pdf",
#       height = 131.25, width = 175, dpi = 300, units = "mm")

#### RNA-seq plants plots ####
# total plant mass was not (could not be) measured
# for plants selected for RNA-seq 

rna_height_plot <- ggplot(phenotype_data %>% filter(RNA_SEQ==1),
                      aes(x=HABITAT, y=HEIGHT,
                          shape=HABITAT, fill=HABITAT)) +
  
  geom_point(size=3, alpha=.75,
             position = position_jitter(width = .2)) +
  geom_violin(alpha=0.25) +
  theme_bw(base_size = 8) +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("gold2", "forestgreen")) +
  scale_color_manual(values = c("gold2", "forestgreen")) +
  scale_shape_manual(values=c(24,21)) +
  labs(x="", y="Height (cm)")

rna_leaf_num_plot <- ggplot(phenotype_data %>% filter(RNA_SEQ==1),
                        aes(x=HABITAT, y=LEAF_num,
                            shape=HABITAT, fill=HABITAT)) +
  geom_point(size=3, alpha=.75,
             position = position_jitter(width = .2, height = 0)) +
  geom_violin(alpha=0.25) +
  theme_bw(base_size = 8) +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("gold2", "forestgreen")) +
  scale_color_manual(values = c("gold2", "forestgreen")) +
  scale_shape_manual(values=c(24,21)) +
  labs(x="Ecotype", y="Num. of leaves")

rna_leaf_mass_dry_plot <- ggplot(phenotype_data %>% filter(RNA_SEQ==1),
                             aes(x=HABITAT, y=LEAF_MASS_DRY,
                                 shape=HABITAT, fill=HABITAT)) +
  
  geom_point(size=3, alpha=.75,
             position = position_jitter(width = .2)) +
  geom_violin(alpha=0.25) +
  theme_bw(base_size = 8) +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("gold2", "forestgreen")) +
  scale_color_manual(values = c("gold2", "forestgreen")) +
  scale_shape_manual(values=c(24,21)) +
  labs(x="Ecotype", y="Leaf mass (g)")

rna_seedling_trait_plot <- (rna_height_plot | rna_leaf_num_plot) / (rna_leaf_mass_dry_plot | plot_spacer())

#### linear models/t-tests ####
# relevel so non-dune is the reference/intercept
phenotype_data$HABITAT <- as.factor(phenotype_data$HABITAT)
phenotype_data$HABITAT <- relevel(phenotype_data$HABITAT, ref = 2)

fit_height <- lm(HEIGHT ~ HABITAT, data = phenotype_data)
summary(fit_height)
t.test(HEIGHT ~ HABITAT, data = phenotype_data)

fit_leaf_num <- lm(LEAF_num ~ HABITAT, data=phenotype_data %>% filter(RNA_SEQ==1))
summary(fit_leaf_num)
t.test(LEAF_num ~ HABITAT, data=phenotype_data %>% filter(RNA_SEQ==1))
mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
}
tmp <- phenotype_data %>% 
    filter(HABITAT=="Non-dune") %>% 
    dplyr::select(LEAF_num) %>% 
    na.omit()
mode(tmp$LEAF_num)


fit_leaf_mass <- lm(LEAF_MASS_DRY ~ HABITAT, data = phenotype_data)
summary(fit_leaf_mass)
t.test(LEAF_MASS_DRY ~ HABITAT, data = phenotype_data)

fit_plant_mass <- lm(PLANT_MASS ~ HABITAT, data = phenotype_data)
summary(fit_plant_mass)
t.test(PLANT_MASS ~ HABITAT, data = phenotype_data)

#fit_leaf_mass_proportion <- lm(LEAF_MASS_DRY/PLANT_MASS ~ HABITAT, data = phenotype_data)
#summary(fit_leaf_mass_proportion)
