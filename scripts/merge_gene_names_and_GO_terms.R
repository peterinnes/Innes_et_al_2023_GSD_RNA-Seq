#### parsing BLAST hits for GO analysis ####
library(dplyr)
library(ggplot2)
library(viridis)
library(ggrepel)
library(patchwork)

# read-in file that has BLASTx matches between Ha412 genes and Araport11 genes
Ha412_Ath_mappings <- read.table("analysis/BLAST_out/HAN412_vs_Araport11_1e-20.txt",
                                 col.names = c("Ha412_gene", "Ath_gene"))
Ha412_Ath_mappings$Ha412_gene <- gsub("mRNA","gene", Ha412_Ath_mappings$Ha412_gene)
#Ha412_Ath_mappings <- read.table("analysis/BLAST_out/top_hits_HAN412_vs_Araport11.txt", col.names = c("Ha412_gene", "Ath_gene"))
#Ha412_Ath_mappings <- read.table("analysis/BLAST_out/temp", col.names = c("Ha412_gene", "Ath_gene"))
#stringtie_Ath_mappings <- read.table("analysis/BLAST_out/top_hits_stringtie_vs_Araport11.txt", col.names = c("stringtie_gene", "Ath_gene"))
#stringtie_Ha412_mappings <- read.table("data/stringtie/merged_transcripts.stringtie.Ha412_gene_names.txt", col.names = c("stringtie_gene", "Ha412_gene"))

# read-in our re-formatted GO associations file 
# (downloaded from arabidopsis.org and reformatted using parse_GO_terms.py)
Ath_GO_associations <- read.table("data/ref_genome_Araport11/ATH_GO_GOSLIM.id2gos.txt",
                                  sep = '\t', col.names = c("Ath_gene", "GO_terms"))
# the other Ath Go associations file... unclear how this differs from the 
# GOSLIM associations
Ath_GAF <- read.table("data/ref_genome_Araport11/gene_association.tair.gaf",
                      skip = 5, sep = '\t')
  
# merge by Ath gene name to get GO terms for Ha412 genes.
Ha412_GO_associations <- left_join(Ha412_Ath_mappings, Ath_GO_associations,
                                   by="Ath_gene") %>%
  dplyr::select(Ha412_gene, GO_terms)

## get rid of the prefix "mRNA" so it's just the gene id.
#Ha412_GO_associations$Ha412_gene <- gsub("mRNA:", "",
#Ha412_GO_associations$Ha412_gene)

head(Ha412_GO_associations)
dim(Ha412_GO_associations)

# write to file
write.table(Ha412_GO_associations,
            file = "analysis/GO_analysis/Ha412_GO_GOSLIM.id2gos.txt",
            quote = F, sep = "\t", row.names = F, col.names = F)

#### visualize GO enrichment results ####
# read-in results from GOMCL clustering
GO_results_DS <- read.table("analysis/GO_analysis/results_DS_rMATS_GO_enrichment_GOsize3000_OC_Ct0.5I1.5.clstr",
                         header = T, sep = "\t") %>%
  mutate(fold_enrichment=(x.cats.test/X.total.test)/(n.cats.ref/N.total.ref)) %>%
  group_by(Clstr) %>%
  slice(which.max(x.cats.test)) #keep only one gene per cluster

## Code for reading output directly from GOATOOLs (i.e. without clustering)
# <- read.table("analysis/GO_analysis/results_DS_rMATS_GO_enrichment.txt",
#                  header = T, sep = '\t') %>%
#  mutate(ratio_in_study=sapply(ratio_in_study, function(x) eval(parse(text=x))),
#         ratio_in_pop=sapply(ratio_in_pop,
#                             function(x) eval(parse(text=x))),
#         FC=ratio_in_study/ratio_in_pop) %>%
#  filter(enrichment=="e", NS="BP")


GO_results_DE_dune <- read.table("analysis/GO_analysis/results_DE_dune_noLFCthreshold_GO_enrichment_GOsize3000_OC_Ct0.5I1.5.clstr",
                                 header = T, sep = '\t') %>%
  mutate(fold_enrichment=(x.cats.test/X.total.test)/(n.cats.ref/N.total.ref)) %>%
  group_by(Clstr) %>%
  dplyr::slice(which.max(x.cats.test)) #keep only one gene per cluster

# what is the degree of overlap betwen the 'response to abscisic acid' and
# 'response to water deprivation'? 50 out out 138 and 132 genes respectively
tmp <- unlist(strsplit(GO_results_DE_dune$Genes.in.test.set[1], split = "\\Q|\\E"))
tmp2 <- unlist(strsplit(GO_results_DE_dune$Genes.in.test.set[3], split = "\\Q|\\E"))
length(intersect(tmp,tmp2))


GO_results_DE_non_dune <- read.table("analysis/GO_analysis/results_DE_non-dune_noLFCthreshold_GO_enrichment_GOsize3000_OC_Ct0.5I1.5.clstr",
                                     header = T, sep = '\t') %>%
  mutate(fold_enrichment=(x.cats.test/X.total.test)/(n.cats.ref/N.total.ref)) %>%
  group_by(Clstr) %>%
  slice(which.max(x.cats.test)) #keep only one gene per cluster
  

# plotting
plot_GO_DS <- ggplot(data=GO_results_DS, aes(x=fold_enrichment,
                                                      y=-log10(adj.p.value))) +
  
  geom_point(aes(size=x.cats.test),
             alpha=0.9, color="#4D4D4D") +
  scale_size(range = c(5,15), breaks = c(20,40,80), name = "Num. genes") +
  
  geom_label_repel(aes(x=fold_enrichment,y=-log10(adj.p.value),
                       label = Description),
                   inherit.aes = T, size = 4.5, label.size = NA,
                   point.padding = 10, label.padding = 0) +
  
  theme_bw() +
  theme(text = element_text(size=24),
        #legend.position = c(.85,.75),
        legend.position = c(.875,.65),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        legend.background = element_blank(),
        legend.box.background = element_rect(color = "black"),
        #panel.grid.minor = element_blank(),
        #panel.grid.major = element_blank(),
        plot.title = element_text(size=18),
        axis.text = element_text(size=18)) +
  labs(x="Fold Enrichment")#, title = "DS")
plot_GO_DS

plot_GO_DE_dune <- ggplot(data=GO_results_DE_dune,
                              aes(x=fold_enrichment,y=-log10(adj.p.value))) +
  
  geom_point(aes(size=x.cats.test),
             alpha=0.9, color="#4D4D4D") +
  #scale_color_viridis() +
  scale_size(range = c(5,20), breaks = c(25, 50, 100), name = "Num. genes") +
  #geom_label(label=GO_results_DE_dune$Description) +
  #geom_label_repel(aes(x=fold_enrichment,y=-log10(adj.p.value),
  #                     label = Description),
  #                 inherit.aes = T, size = 4.5, label.size = NA,
  #                 point.padding = 10, label.padding = 0) +
  
  theme_bw() +
  theme(text = element_text(size=24),
        #legend.position = c(.85,.75),
        legend.position = c(.875,.65),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        legend.background = element_blank(),
        legend.box.background = element_rect(color = "black"),
        #panel.grid.minor = element_blank(),
        #panel.grid.major = element_blank(),
        plot.title = element_text(size=18),
        axis.text = element_text(size=18)) +
  scale_y_continuous(breaks = seq(0, 6, by = 1)) +
  labs(x="Fold Enrichment")#, title = "DE (dune up-regulated)")
plot_GO_DE_dune

plot_GO_DE_non_dune <- ggplot(data=GO_results_DE_non_dune,
                              aes(x=fold_enrichment,
                                        y=-log10(adj.p.value))) +
  geom_point(aes(size=x.cats.test),
             alpha=0.9, color="#4D4D4D") +
  #geom_text(aes(label=ifelse(x.cats.test>40, Description,''))) +
  scale_size(range = c(5,20), breaks = c(25,50,100), name = "Num. genes") +
  geom_label_repel(aes(x=fold_enrichment,y=-log10(adj.p.value), 
                       label=ifelse(x.cats.test>20, Description,'')),
                   max.overlaps = 12, inherit.aes = T, size = 4.5,
                   label.size = NA, point.padding = 10, label.padding = 0) +
  theme_bw() +
  theme(text = element_text(size=24),
        legend.position = c(.9,.65),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        legend.background = element_blank(),
        legend.box.background = element_rect(color = "black"),
        #panel.grid.minor = element_blank(),
        #panel.grid.major = element_blank(),
        plot.title = element_text(size=18),
        axis.text = element_text(size=18)) +
  labs(x="Fold Enrichment")#, title = "DE (non-dune up-regulated)")
plot_GO_DE_non_dune

GO_enrich_fig <- plot_GO_DE_dune / plot_GO_DS
ggsave(filename = "figures/GO_enrich_fig.png", device = "png",
       plot=GO_enrich_fig, height=9, width=7.5, dpi = 300,
       units = "in", bg = "white")

ggsave(filename = "figures/plot_GO_DS.png", device = "png",
       plot=plot_GO_DS, height=4.5, width=6, dpi = 300,
       units = "in", bg = "white")

ggsave(filename = "figures/plot_GO_DE_dune.png", device = "png",
       plot=plot_GO_DE_dune, height=4.5, width=6, dpi = 300,
       units = "in", bg = "white")

#### examine DE–DS intersect genes ####
GO_results_intersect_DE_DS <- read.table("analysis/GO_analysis/results_intersect_DE–DS_noLFCthreshold_GO_enrichment.txt", header=T,
                                         sep="\t")
# list the DE–DS "embryo development..." genes. Notably, ZERO
# of them are on chrom 11. 
intersect_DE_DS_embryo_dev_genes <- 
  unlist(strsplit(GO_results_intersect_DE_DS$study_items[1], split = ", "))

# list the DS "embryo development..." genes
DS_embryo_dev_genes <- read.table("analysis/GO_analysis/results_DS_rMATS_GO_enrichment.txt",
                                  header=T, sep="\t") %>%
  filter(name=="embryo development ending in seed dormancy") %>%
  dplyr::select(study_items) 

DS_embryo_dev_genes <- unlist(strsplit(x = DS_embryo_dev_genes$study_items,split = ", "))

# read-in GO term clustering results
GOMCL_results <- read.table("analysis/GO_analysis/results_intersect_DE–DS_noLFCthreshold_GO_enrichment_GOsize3000_OC_Ct0.5I1.5.BP_CC_MF.clstr", header=T, sep="\t") %>%
  mutate(fold_enrichment=(x.cats.test/X.total.test)/(n.cats.ref/N.total.ref)) %>%
  group_by(Type, Clstr) %>%
  slice(which.max(n.cats.ref)) %>%
  ungroup()

plot_GO_intersect_DE_DS <- ggplot(data=GOMCL_results, aes(x=fold_enrichment,
                                             y=-log10(adj.p.value))) +
  
  geom_point(aes(size=x.cats.test, color=Type),
             alpha=0.9) +
  #scale_color_viridis() +
  scale_size(range = c(5,15), breaks = c(5,10,30), name = "Num. genes") +
  scale_color_grey(name = "GO type",
                    labels = c("Biological process",
                                                 "Cellular component",
                                                 "Molecular function")) +
  
  geom_label_repel(aes(x=fold_enrichment,y=-log10(adj.p.value),
                       label = Description),
                   inherit.aes = T, size = 4.5, label.size = NA, 
                   point.padding = 10, label.padding = 0) +
  
  theme_bw() +
  theme(text = element_text(size=24),
        legend.position = c(.8,.5),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        legend.background = element_blank(),
        legend.box.background = element_rect(color = "black"),
        #panel.grid.minor = element_blank(),
        #panel.grid.major = element_blank(),
        plot.title = element_text(size=18),
        axis.text = element_text(size=18)) +
  guides(color = guide_legend(override.aes = list(size=5))) +
  labs(x="Fold Enrichment")#, title = "DE–DS")

plot_GO_intersect_DE_DS
ggsave(filename = "figures/plot_GO_intersect_DE_DS.png", device = "png",
       plot=plot_GO_intersect_DE_DS, height=4.5, width=6, dpi = 300,
       units = "in", bg = "white")

#### Graveyard ####
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