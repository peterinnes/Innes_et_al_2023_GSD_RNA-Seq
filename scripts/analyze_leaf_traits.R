#### read-in phenotype data ####
a <- read.csv("data2/GSDgreenhouseExpt2018_Data_sheet1.csv", header = T) %>%
  dplyr::select(HABITAT, POPULATION, plant_number, GERM_DATE, ACHENE_LEN, HYPO_LEN, ROOT_LEN,
                HEIGHT, LEAF_num, LEAF_MASS_WET, LEAF_MASS_DRY, RNA_EXT)
b <- read.csv("data2/GSDgreenhouseExpt2018_Data_sheet2.csv", header = T) %>%
  dplyr::select(HABITAT, POPULATION, plant_number, code_number, RNA_SEQ)

phenotype_data <- inner_join(a,b) %>%
  filter(RNA_SEQ==1)


ggplot(data=phenotype_data, aes(x=HABITAT, y=LEAF_MASS_DRY)) +
  geom_violin()

fit_leaf_mass <- lm(LEAF_MASS_DRY ~ -1 + HABITAT, data = phenotype_data)
summary(fit_leaf_mass)