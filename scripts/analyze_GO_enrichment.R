#### analyze GO enrichment results ####
library(forcats)

#### Dot plot visulations ####

#### DE dune-upregulated genes ####
GO_results_DE_dune_BP <- read.table("analysis/GO_analysis/GOMCL_BP/results_DE_dune_noLFCthreshold_GO_enrichment_GOsize3000_OC_Ct0.5I1.5.clstr",
                                                          header = T, sep = '\t') %>%
  mutate(fold_enrichment=(x.cats.test/X.total.test)/(n.cats.ref/N.total.ref)) %>%
  group_by(Clstr) %>%
  dplyr::slice(which.max(x.cats.test)) #keep only one term per cluster, the one with greatest num genes

# only two CC GO terms...somehow only one made it into the clustering results.
GO_results_DE_dune_CC <- read.table("analysis/GO_analysis/GOMCL_CC/results_DE_dune_noLFCthreshold_GO_enrichment_GOsize3000_OC_Ct0.5I1.5.clstr",
                                    header = T, sep = '\t') %>%
  mutate(fold_enrichment=(x.cats.test/X.total.test)/(n.cats.ref/N.total.ref)) %>%
  group_by(Clstr) %>%
  dplyr::slice(which.max(x.cats.test))
GO_results_DE_dune_CC$Description[1] <- "...endoplasmic reticulum membrane"

GO_results_DE_dune_MF <- read.table("analysis/GO_analysis/GOMCL_MF/results_DE_dune_noLFCthreshold_GO_enrichment_GOsize3000_OC_Ct0.5I1.5.clstr",
                                    header = T, sep = '\t') %>%
  mutate(fold_enrichment=(x.cats.test/X.total.test)/(n.cats.ref/N.total.ref)) %>%
  group_by(Clstr) %>%
  dplyr::slice(which.max(x.cats.test))

GO_results_DE_dune_ALL <- bind_rows(GO_results_DE_dune_BP, GO_results_DE_dune_CC) %>%
  bind_rows(GO_results_DE_dune_MF) %>%
  mutate(Type_num=gsub("BP", 1, Type)) %>%
  mutate(Type_num=gsub("CC", 2, Type_num)) %>%
  mutate(Type_num=gsub("MF", 3, Type_num)) %>%
  mutate(Set="DE dune upregulated")

GO_DE_dune_dotplot <- ggplot(data = GO_results_DE_dune_ALL,
       aes(x=fold_enrichment, y=fct_reorder(Description, -log10(p.value),
                                            .desc=F))) +
  geom_point(aes(size=-log10(p.value), color=Type)) +
  scale_color_manual(name = "GO category", values = c("black", "grey50", "grey75"),
                     labels = c("Biological process",
                                "Cellular component",
                                "Molecular function"), guide = "none") +
  scale_size_continuous(name = bquote(paste("-",log[10], "(p-value)")),
                        range = c(2,6), breaks = c(2,4,8,16)) + 
  
  #scale_x_continuous(expand = c(.2,0)) +
  #coord_fixed(ratio = 1.75) +
  #scale_x_continuous(limits = c(1,9)) +
  #scale_y_discrete(expand = c(0.4,.4)) +
  coord_fixed(1.25) +
  theme_bw(base_size = 18) +
  #theme(axis.text = element_text(size=))
  theme(legend.position = "none") + #c(.8,.8),
        #text=element_text(size=18),
        #axis.text = element_text(size=12),
        #legend.text = element_text(size=12),
        #legend.title = element_text(size=12)) +
        #legend.background = element_blank(),
        #legend.box.background = element_rect(color = "black"),
        #plot.title = element_text(size=18)) +
  guides(color = guide_legend(override.aes = list(size=4))) +
  labs(y="", x="Fold enrichment")

GO_DE_dune_dotplot

png(filename = "figures/GO_DE_dune_dotplot.png",width = 6, height = 4.5,
    units = "in", res = 300)
GO_DE_dune_dotplot
dev.off()

#### DE non-dune upregulated genes ####
GO_results_DE_non_dune_BP <- read.table("analysis/GO_analysis/GOMCL_BP/results_DE_non-dune_noLFCthreshold_GO_enrichment_GOsize3000_OC_Ct0.5I1.5.clstr",
                                    header = T, sep = '\t') %>%
  mutate(fold_enrichment=(x.cats.test/X.total.test)/(n.cats.ref/N.total.ref)) %>%
  group_by(Clstr) %>%
  dplyr::slice(which.max(x.cats.test)) #keep only one term per cluster, the one with greatest num genes

# only two CC GO terms...somehow only one made it into the clustering results.
GO_results_DE_non_dune_CC <- read.table("analysis/GO_analysis/GOMCL_CC/results_DE_non-dune_noLFCthreshold_GO_enrichment_GOsize3000_OC_Ct0.5I1.5.clstr",
                                    header = T, sep = '\t') %>%
  mutate(fold_enrichment=(x.cats.test/X.total.test)/(n.cats.ref/N.total.ref)) %>%
  group_by(Clstr) %>%
  dplyr::slice(which.max(x.cats.test))

GO_results_DE_non_dune_MF <- read.table("analysis/GO_analysis/GOMCL_MF/results_DE_non-dune_noLFCthreshold_GO_enrichment_GOsize3000_OC_Ct0.5I1.5.clstr",
                                    header = T, sep = '\t') %>%
  mutate(fold_enrichment=(x.cats.test/X.total.test)/(n.cats.ref/N.total.ref)) %>%
  group_by(Clstr) %>%
  dplyr::slice(which.max(x.cats.test))

GO_results_DE_non_dune_ALL <- bind_rows(GO_results_DE_non_dune_BP,
                                        GO_results_DE_non_dune_CC) %>%
  bind_rows(GO_results_DE_non_dune_MF) %>%
  mutate(Type_num=gsub("BP", 1, Type)) %>%
  mutate(Type_num=gsub("CC", 2, Type_num)) %>%
  mutate(Type_num=gsub("MF", 3, Type_num)) %>%
  mutate(Set="DE non-dune upregulated")

# what is the degree of overlap between the 'embryo dev. ending in seed dorm.' and
# 'seed development'? only 3 out out 135 and 70 genes respectively
tmp <- unlist(strsplit(GO_results_DE_non_dune_ALL$Genes.in.test.set[2], split = "\\Q|\\E"))
tmp2 <- unlist(strsplit(GO_results_DE_non_dune_ALL$Genes.in.test.set[5], split = "\\Q|\\E"))
length(intersect(tmp,tmp2))

# dot plot
GO_DE_non_dune_dotplot <- ggplot(data = GO_results_DE_non_dune_ALL,
                             aes(x=fold_enrichment, y=fct_reorder(Description, -log10(p.value),
                                                                  .desc=F))) +
  geom_point(aes(size=-log10(p.value), color=Type)) +
  scale_color_manual(name = "GO category", values = c("black", "grey50", "grey75"),
                     labels = c("Biological process",
                                "Cellular component",
                                "Molecular function"), guide = "none") +
  scale_size_continuous(name = bquote(paste("-",log[10], "(p-value)")),
                        range = c(2,10), breaks = c(2,4,16,72)) + 
  
  theme_bw(base_size = 12) +
  theme(axis.text.y = element_text(size=8),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8)) +
  guides(color = guide_legend(override.aes = list(size=3))) +
  labs(y="", x="Fold enrichment")

ggsave(filename = "figures/GO_DE_non_dune_dotplot.pdf",
       plot = GO_DE_non_dune_dotplot, device = "pdf", width = 175, height = 232.75,
       units = "mm")

png(filename = "figures/GO_DE_non_dune_dotplot.png",width = 10, height = 15,
    units = "in", res = 300)
GO_DE_non_dune_dotplot
dev.off()

#### DS genes #### 
GO_results_DS_dune_BP <- read.table("analysis/GO_analysis/GOMCL_BP/results_DS_rMATS_GO_enrichment_GOsize3000_OC_Ct0.5I1.5.clstr",
                                    header = T, sep = '\t') %>%
  mutate(fold_enrichment=(x.cats.test/X.total.test)/(n.cats.ref/N.total.ref)) %>%
  group_by(Clstr) %>%
  dplyr::slice(which.max(x.cats.test)) #keep only one term per cluster, the one with greatest num genes

# only two CC GO terms...somehow only one made it into the clustering results.
GO_results_DS_dune_CC <- read.table("analysis/GO_analysis/GOMCL_CC/results_DS_rMATS_GO_enrichment_GOsize3000_OC_Ct0.5I1.5.clstr",
                                    header = T, sep = '\t') %>%
  mutate(fold_enrichment=(x.cats.test/X.total.test)/(n.cats.ref/N.total.ref)) %>%
  group_by(Clstr) %>%
  #dplyr::slice(which.min(adj.p.value))
  dplyr::slice(which.max(x.cats.test))

GO_results_DS_dune_MF <- read.table("analysis/GO_analysis/GOMCL_MF/results_DS_rMATS_GO_enrichment_GOsize3000_OC_Ct0.5I1.5.clstr",
                                    header = T, sep = '\t') %>%
  mutate(fold_enrichment=(x.cats.test/X.total.test)/(n.cats.ref/N.total.ref)) %>%
  group_by(Clstr) %>%
  dplyr::slice(which.max(x.cats.test))

GO_results_DS_ALL <- bind_rows(GO_results_DS_dune_BP, GO_results_DS_dune_CC) %>%
  bind_rows(GO_results_DS_dune_MF) %>%
  mutate(Type_num=gsub("BP", 1, Type)) %>%
  mutate(Type_num=gsub("CC", 2, Type_num)) %>%
  mutate(Type_num=gsub("MF", 3, Type_num)) %>%
  mutate(Set="DS")

GO_results_DS_ALL$Description[5] <- "glycine decarboxylation..."
GO_results_DS_ALL$Description[2] <- "embryo dev. ending in seed dorm."

GO_DS_dotplot <- ggplot(data = GO_results_DS_ALL,
                        aes(x=fold_enrichment,
                            y=fct_reorder(Description, -log10(adj.p.value),
                                          .desc=F))) +
  
  geom_point(aes(size=-log10(p.value), color=Type)) +
  scale_color_manual(name = "GO category", values = c("black", "grey50", "grey75"),
                     labels = c("Biological process",
                                "Cellular component",
                                "Molecular function")) +
  scale_size_continuous(name = bquote(paste("-", log[10], "(p-value)")),
                        range = c(2,10), breaks = c(2,4,8,16)) + 
  
  scale_x_continuous(expand = c(.2,0)) +
  coord_fixed(2) +
  
  theme_bw(base_size=18) +
  theme(axis.text.y = element_text(size=8),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8)) +
  guides(color = guide_legend(override.aes = list(size=3))) +
  labs(x="Fold enrichment", y="")

GO_DS_dotplot

png(filename = "figures/GO_DS_dotplot.png", res = 300, width = 8,
    height = 4.5, units = "in")
GO_DS_dotplot
dev.off()

#### DE–DS intersect genes ####
GO_results_intersect_BP <- read.table("analysis/GO_analysis/GOMCL_BP/results_intersect_DE–DS_noLFCthreshold_GO_enrichment_GOsize3000_OC_Ct0.5I1.5.clstr",
                                      header=T, sep="\t") %>%
  mutate(fold_enrichment=(x.cats.test/X.total.test)/(n.cats.ref/N.total.ref)) %>%
  group_by(Clstr) %>%
  dplyr::slice(which.max(x.cats.test))

GO_results_intersect_CC <- read.table("analysis/GO_analysis/GOMCL_CC/results_intersect_DE–DS_noLFCthreshold_GO_enrichment_GOsize3000_OC_Ct0.5I1.5.clstr",
                                      header=T, sep="\t") %>%
  mutate(fold_enrichment=(x.cats.test/X.total.test)/(n.cats.ref/N.total.ref)) %>%
  group_by(Clstr) %>%
  dplyr::slice(which.max(x.cats.test))

GO_results_intersect_MF <- read.table("analysis/GO_analysis/GOMCL_MF/results_intersect_DE–DS_noLFCthreshold_GO_enrichment_GOsize3000_OC_Ct0.5I1.5.clstr",
                                      header=T, sep="\t") %>%
  mutate(fold_enrichment=(x.cats.test/X.total.test)/(n.cats.ref/N.total.ref)) %>%
  group_by(Clstr) %>%
  dplyr::slice(which.max(x.cats.test))

GO_results_intersect_ALL <- bind_rows(GO_results_intersect_BP,
                                      GO_results_intersect_CC) %>%
  bind_rows(GO_results_intersect_MF) %>%
  mutate(Type_num=gsub("BP", 1, Type)) %>%
  mutate(Type_num=gsub("CC", 2, Type_num)) %>%
  mutate(Type_num=gsub("MF", 3, Type_num)) %>%
  mutate(Set="intersect")

GO_results_intersect_ALL$Description[1] <- "embryo dev. ending in seed dorm."

save(GO_results_DE_dune_ALL, GO_results_DS_ALL, GO_results_intersect_ALL,
     file = "data2/Rdata/GO_results_main_text.Rdata")

GO_intersect_dotplot <- ggplot(data = GO_results_intersect_ALL,
                        aes(x=fold_enrichment,
                            y=fct_reorder(Description, -log10(p.value),
                                          .desc=F))) +
  
  geom_point(aes(size=-log10(p.value), color=Type)) +
  scale_color_manual(name = "GO category", values = c("black", "grey50", "grey75"),
                     labels = c("Biological process",
                                "Cellular component",
                                "Molecular function")) +
  scale_size_continuous(name = bquote(paste("-", log[10], "(p-value)")),
                        range = c(2,8), breaks = c(2,4,8,16)) + 
  
  scale_x_continuous(expand = c(.1,0)) +
  scale_y_discrete(expand = c(.25,0), position = "right") +
  coord_fixed(1) +
  
  theme_bw(base_size = 18) +
  theme(#legend.position = "left", 
        #legend.text = element_text(size=12)) +
        legend.position = "none") +
        #legend.title = element_text(size=8)) +
  guides(color = guide_legend(override.aes = list(size=4))) +
  labs(x="Fold enrichment", y="")

png(filename = "figures/GO_intersect_dotplot.png", res = 300, width = 6,
    height = 4.5, units = "in")
GO_intersect_dotplot
dev.off()

# plot DE dune, DS, intersect sets together with facet_wrap()
tmp <- rbind(GO_results_DE_dune_ALL, GO_results_DS_ALL, GO_results_intersect_ALL)
GO_dotplot_facet <- ggplot(data = tmp,
                               aes(x=fold_enrichment,
                                   y=fct_reorder(Description, -log10(p.value),
                                                 .desc=F))) +
  facet_wrap(~Set) +
  geom_point(aes(size=-log10(p.value), color=Type)) +
  scale_color_manual(name = "GO category", values = c("black", "grey50", "grey75"),
                     labels = c("Biological process",
                                "Cellular component",
                                "Molecular function")) +
  scale_size_continuous(name = bquote(paste("-", log[10], "(p-value)")),
                        range = c(2,16), breaks = c(2,4,8,16)) + 
  
  #scale_x_continuous(expand = c(.1,0)) +
  #scale_y_discrete(expand = c(.25,0), position = "right") +
  #coord_fixed(1) +
  
  theme_bw(base_size = 14) +
  theme(#legend.position = "left", 
    #legend.text = element_text(size=12)) +
    legend.position = "top") +
  #legend.title = element_text(size=8)) +
  guides(color = guide_legend(override.aes = list(size=4))) +
  labs(x="Fold enrichment", y="")

GO_dotplot_facet

#### look at genes of specific GO terms ####

# list the DE–DS "embryo development..." genes. Notably, ZERO
# of them are on chrom 11. 
intersect_embryo_dev_genes <- 
  unlist(strsplit(GO_results_intersect_ALL$study_items[1], split = ", "))

# list the DE (non-dune upregulated) "embryo development..." and "seed development" genes
DE_embryo_seed_dev_genes <- read.table("analysis/GO_analysis/results_DE_non-dune_noLFCthreshold_GO_enrichment.txt", header=T, sep="\t") %>%
  filter(name %in% c("seed development", "embryo development ending in seed dormancy")) %>%
  dplyr::select(study_items) 
DE_embryo_seed_dev_genes <- unlist(strsplit(x = DE_embryo_seed_dev_genes$study_items,split = ", "))

# list the DS "embryo development..." genes
DS_embryo_dev_genes <- read.table("analysis/GO_analysis/results_DS_rMATS_GO_enrichment.txt",
                                  header=T, sep="\t") %>%
  filter(name=="embryo development ending in seed dormancy") %>%
  dplyr::select(study_items) 
DS_embryo_dev_genes <- unlist(strsplit(x = DS_embryo_dev_genes$study_items,split = ", "))


# list the DS genes related to nitrogen assimilation
DS_nitrogen_assim_genes <- read.table("analysis/GO_analysis/results_DS_rMATS_GO_enrichment.txt",
                                      header=T, sep="\t") %>%
  filter(name=="ammonia assimilation cycle") %>%
  dplyr::select(study_items)
DS_nitrogen_assim_genes <- unlist(strsplit(x = DS_nitrogen_assim_genes$study_items,split = ", "))


#### Graveyard ####
# read-in results from GOMCL clustering
#GO_results_DS <- read.table("analysis/GO_analysis/results_DS_rMATS_GO_enrichment_GOsize3000_OC_Ct0.5I1#.5.clstr",
#                            header = T, sep = "\t") %>%
#  mutate(fold_enrichment=(x.cats.test/X.total.test)/(n.cats.ref/N.total.ref)) %>%
#  group_by(Clstr) %>%
#  slice(which.max(x.cats.test)) #keep only one term per cluster

## Code for reading output directly from GOATOOLs (i.e. without clustering)
# <- read.table("analysis/GO_analysis/results_DS_rMATS_GO_enrichment.txt",
#                  header = T, sep = '\t') %>%
#  mutate(ratio_in_study=sapply(ratio_in_study, function(x) eval(parse(text=x))),
#         ratio_in_pop=sapply(ratio_in_pop,
#                             function(x) eval(parse(text=x))),
#         FC=ratio_in_study/ratio_in_pop) %>%
#  filter(enrichment=="e", NS="BP")


#GO_results_DE_dune <- read.table("analysis/GO_analysis/GOMCL_BP#/results_DE_dune_noLFCthreshold_GO_enrichment_GOsize3000_OC_Ct0.5I1.5.clstr",
#                                 header = T, sep = '\t') %>%
#  mutate(fold_enrichment=(x.cats.test/X.total.test)/(n.cats.ref/N.total.ref)) %>%
#  group_by(Clstr) %>%
#  dplyr::slice(which.max(x.cats.test)) #keep only one gene per cluster
#
## what is the degree of overlap betwen the 'response to abscisic acid' and
## 'response to water deprivation'? 50 out out 138 and 132 genes respectively
#tmp <- unlist(strsplit(GO_results_DE_dune$Genes.in.test.set[1], split = "\\Q|\\E"))
#tmp2 <- unlist(strsplit(GO_results_DE_dune$Genes.in.test.set[3], split = "\\Q|\\E"))
#length(intersect(tmp,tmp2))
#
#
#GO_results_DE_non_dune <- read.table("analysis/GO_analysis/results_DE_non#-dune_noLFCthreshold_GO_enrichment_GOsize3000_OC_Ct0.5I1.5.clstr",
#                                     header = T, sep = '\t') %>%
#  mutate(fold_enrichment=(x.cats.test/X.total.test)/(n.cats.ref/N.total.ref)) %>%
#  group_by(Clstr) %>%
#  slice(which.max(x.cats.test)) #keep only one gene per cluster
#
#
## plotting
#plot_GO_DS <- ggplot(data=GO_results_DS, aes(x=fold_enrichment,
#                                             y=-log10(adj.p.value))) +
#  
#  geom_point(aes(size=x.cats.test),
#             alpha=0.9, color="#4D4D4D") +
#  scale_size(range = c(5,15), breaks = c(20,40,80), name = "Num. genes") +
#  
#  geom_label_repel(aes(x=fold_enrichment,y=-log10(adj.p.value),
#                       label = Description),
#                   inherit.aes = T, size = 4.5, label.size = NA,
#                   point.padding = 10, label.padding = 0) +
#  
#  theme_bw() +
#  theme(text = element_text(size=24),
#        #legend.position = c(.85,.75),
#        legend.position = c(.875,.65),
#        legend.text = element_text(size=12),
#        legend.title = element_text(size=12),
#        legend.background = element_blank(),
#        legend.box.background = element_rect(color = "black"),
#        #panel.grid.minor = element_blank(),
#        #panel.grid.major = element_blank(),
#        plot.title = element_text(size=18),
#        axis.text = element_text(size=18)) +
#  labs(x="Fold Enrichment")#, title = "DS")
#plot_GO_DS
#
#plot_GO_DE_dune <- ggplot(data=GO_results_DE_dune,
#                          aes(x=fold_enrichment,y=-log10(adj.p.value))) +
#  
#  geom_point(aes(size=x.cats.test),
#             alpha=0.9, color="#4D4D4D") +
#  #scale_color_viridis() +
#  scale_size(range = c(5,20), breaks = c(25, 50, 100), name = "Num. genes") +
#  #geom_label(label=GO_results_DE_dune$Description) +
#  #geom_label_repel(aes(x=fold_enrichment,y=-log10(adj.p.value),
#  #                     label = Description),
#  #                 inherit.aes = T, size = 4.5, label.size = NA,
#  #                 point.padding = 10, label.padding = 0) +
#  
#  theme_bw() +
#  theme(text = element_text(size=24),
#        #legend.position = c(.85,.75),
#        legend.position = c(.875,.65),
#        legend.text = element_text(size=12),
#        legend.title = element_text(size=12),
#        legend.background = element_blank(),
#        legend.box.background = element_rect(color = "black"),
#        #panel.grid.minor = element_blank(),
#        #panel.grid.major = element_blank(),
#        plot.title = element_text(size=18),
#        axis.text = element_text(size=18)) +
#  scale_y_continuous(breaks = seq(0, 6, by = 1)) +
#  labs(x="Fold Enrichment")#, title = "DE (dune up-regulated)")
#plot_GO_DE_dune
#
#plot_GO_DE_non_dune <- ggplot(data=GO_results_DE_non_dune,
#                              aes(x=fold_enrichment,
#                                  y=-log10(adj.p.value))) +
#  geom_point(aes(size=x.cats.test),
#             alpha=0.9, color="#4D4D4D") +
#  #geom_text(aes(label=ifelse(x.cats.test>40, Description,''))) +
#  scale_size(range = c(5,20), breaks = c(25,50,100), name = "Num. genes") +
#  geom_label_repel(aes(x=fold_enrichment,y=-log10(adj.p.value), 
#                       label=ifelse(x.cats.test>20, Description,'')),
#                   max.overlaps = 12, inherit.aes = T, size = 4.5,
#                   label.size = NA, point.padding = 10, label.padding = 0) +
#  theme_bw() +
#  theme(text = element_text(size=24),
#        legend.position = c(.9,.65),
#        legend.text = element_text(size=12),
#        legend.title = element_text(size=12),
#        legend.background = element_blank(),
#        legend.box.background = element_rect(color = "black"),
#        #panel.grid.minor = element_blank(),
#        #panel.grid.major = element_blank(),
#        plot.title = element_text(size=18),
#        axis.text = element_text(size=18)) +
#  labs(x="Fold Enrichment")#, title = "DE (non-dune up-regulated)")
#plot_GO_DE_non_dune
#
#GO_enrich_fig <- plot_GO_DE_dune / plot_GO_DS
#ggsave(filename = "figures/GO_enrich_fig.png", device = "png",
#       plot=GO_enrich_fig, height=9, width=7.5, dpi = 300,
#       units = "in", bg = "white")
#
#ggsave(filename = "figures/plot_GO_DS.png", device = "png",
#       plot=plot_GO_DS, height=4.5, width=6, dpi = 300,
#       units = "in", bg = "white")
#
#ggsave(filename = "figures/plot_GO_DE_dune.png", device = "png",
#       plot=plot_GO_DE_dune, height=4.5, width=6, dpi = 300,
#       units = "in", bg = "white")

## read-in GO term clustering results
#GOMCL_results <- read.table("analysis/GO_analysis#/results_intersect_DE–DS_noLFCthreshold_GO_enrichment_GOsize3000_OC_Ct0.5I1.5.BP_CC_MF.clstr", header#=T, sep="\t") %>%
#  mutate(fold_enrichment=(x.cats.test/X.total.test)/(n.cats.ref/N.total.ref)) %>%
#  group_by(Type, Clstr) %>%
#  slice(which.max(n.cats.ref)) %>%
#  ungroup()
#
#plot_GO_intersect_DE_DS <- ggplot(data=GOMCL_results, aes(x=fold_enrichment,
#                                                          y=-log10(p.value))) +
#  
#  geom_point(aes(size=x.cats.test, color=Type),
#             alpha=0.9) +
#  #scale_color_viridis() +
#  scale_size(range = c(5,15), breaks = c(5,10,30), name = "Num. genes") +
#  scale_color_grey(name = "GO type",
#                   labels = c("Biological process",
#                              "Cellular component",
#                              "Molecular function")) +
#  
#  geom_label_repel(aes(x=fold_enrichment,y=-log10(p.value),
#                       label = Description),
#                   inherit.aes = T, size = 4.5, label.size = NA, 
#                   point.padding = 10, label.padding = 0) +
#  
#  theme_bw() +
#  theme(text = element_text(size=24),
#        legend.position = c(.8,.5),
#        legend.text = element_text(size=12),
#        legend.title = element_text(size=12),
#        legend.background = element_blank(),
#        legend.box.background = element_rect(color = "black"),
#        #panel.grid.minor = element_blank(),
#        #panel.grid.major = element_blank(),
#        plot.title = element_text(size=18),
#        axis.text = element_text(size=18)) +
#  guides(color = guide_legend(override.aes = list(size=5))) +
#  labs(x="Fold Enrichment")#, title = "DE–DS")
#
#plot_GO_intersect_DE_DS
#ggsave(filename = "figures/plot_GO_intersect_DE_DS.png", device = "png",
#       plot=plot_GO_intersect_DE_DS, height=4.5, width=6, dpi = 300,
#       units = "in", bg = "white")


## Revigo results table of DE genes upregulated in Dune ecotype
#revigo_DE_dune <- read.table("analysis/GO_analysis/Revigo_BP_table.DE_dune.tsv", header = T)
#
## Revigo results table of all DS genes (all unique rMATS + DEXSeq significant DS/DEU genes)
#revigo_DS_rMATS<- read.table("analysis/GO_analysis/Revigo_BP_Table.tsv", header=T, na.strings = "null")
#
#
#plot_revigo_DE_dune <-  ggplot(data=revigo_DE_dune, aes(x = PC_1, y = PC_0, size=Value, color=Value)) +
#  scale_size_continuous(range = c(3,11)) +
#  geom_point(alpha=.5) +
#  geom_text_repel(label=ifelse(revigo_DE_dune$Value>10,as.character(revigo_DE_dune$Name),''), force = 5#, force_pull = 0.5, direction = "both") +
#  labs(x="MDS2 (semantic space x)", y="MDS1 (semantic space y)",
#       title="DE (dune)") +
#  theme_bw() +
#  geom_hline(yintercept = 0, linetype="dashed", color="gray50", alpha=.5) +
#  geom_vline(xintercept = 0, linetype="dashed", color="gray50", alpha=.5) +
#  theme(text = element_text(size=36),
#        legend.position = "none",
#        legend.title = element_blank())
#        #legend.position = "bottom",
#        #legend.margin=margin(0,0,0,0),
#        #legend.box.margin=margin(0,0,0,0))
#plot_revigo_DE_dune
#ggsave(filename = "figures/revigo_mds_DE_dune.png", device = "png", plot=plot_revigo_DE_dune, height=6, #width=8, dpi = 600, units = "in", bg = "white")
#
#plot_revigo_DS <-  ggplot(data=revigo_DS_rMATS, aes(x = PC_1, y = PC_0, size=-log(Value), color=-log#(Value))) +
#  scale_size_continuous(range = c(1,5)) +
#  geom_point(alpha=.5) +
#  geom_text_repel(label=ifelse(revigo_DS_rMATS$Value<1e-2,as.character(revigo_DS_rMATS$Name),''), force #= 5, force_pull = 0.5, direction = "both") +
#  labs(x="MDS2 (semantic space x)", y="",
#       title="DS") +
#  theme_bw() +
#  geom_hline(yintercept = 0, linetype="dashed", color="gray50", alpha=.5) +
#  geom_vline(xintercept = 0, linetype="dashed", color="gray50", alpha=.5) +
#  theme(text = element_text(size=36),
#        legend.position = "none",
#        legend.title = element_blank())
#plot_revigo_DS
##legend.position = "bottom",
##legend.margin=margin(0,0,0,0),
##legend.box.margin=margin(0,0,0,0))
#ggsave(filename = "figures/revigo_mds_DS.png", device = "png", plot=plot_revigo_DS, height=6, width=8, #dpi = 600, units = "in", bg = "white")