library(ggplot2)
library(patchwork)
library(eulerr)
library(ggplotify)

#### Fig 1 (PCA) ####
load("data2/Rdata/gsd_sampling_map.Rdata") #from gsd_sampling_map.R
load("data2/Rdata/snp_pca.Rdata") #from pca.R
load("data2/Rdata/expr_pca.Rdata") #from pca.R
load("data2/Rdata/splice_pca.Rdata") #from pca.R

# Fig 1A
gsd_sampling_map <- ggmap(map) +
  geom_point(data = coordinates,
             mapping = aes(x=long, y=lat, shape=Ecotype, fill=Ecotype),
             size=4, alpha=.75) +
  scale_shape_manual(values = c(24,21)) +
  scale_fill_manual(values=c("gold2", "forestgreen")) +
  theme_bw(base_size = 12) +
  theme(legend.key = element_rect(fill=NA),
        legend.position = "bottom",
        legend.margin=margin(t=-10),
        legend.spacing.x = unit(.1, 'cm'),
        legend.background = element_blank()) +
  scale_x_continuous(expand=c(0,0), breaks = c(-105.6, -105.5)) +
  scale_y_continuous(expand=c(0,0), breaks = c(37.65, 37.7, 37.75)) +
  labs(x="Lon", y="Lat") +
  ggsn::scalebar(data = scale_bar_coords, location = "bottomleft", transform = T,
                 dist = 2, dist_unit = "km", height = .02,
                 st.dist = .075,st.bottom = F, st.size = 3, box.fill = c("black","white"))

ggsave("figures/gsd_sampling_map.pdf", plot = gsd_sampling_map, device = "pdf",
       height = 65.625, width = 87.5, units = "mm", dpi = 300)

# Fig 1B
hull_snp_d <- filter(snp_pca_df, ecotype=="dune") %>%
  dplyr::slice(chull(EV1,EV2))
hull_snp_nd <- filter(snp_pca_df, ecotype=="non-dune") %>%
  dplyr::slice(chull(EV1,EV2))
plot_1B <- ggplot(snp_pca_df,
                       aes(x=EV1, y=EV2, shape=ecotype, fill=ecotype)) +
  # plot each plant and plot convex hulls for each ecotype
  geom_point(size = 4, alpha=.75) +
  geom_polygon(data=hull_snp_d, alpha=.5) +
  geom_polygon(data=hull_snp_nd, alpha=.5) +
  
  labs(x=paste0("PC1 ","(",round(snp_pve_df[1,2], 2),"%)"),
       y=paste0("PC2 ", "(",round(snp_pve_df[2,2],2),"%)")) +
  
  theme_bw(base_size  = 12) +
  theme(legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +

  geom_hline(yintercept = 0, lty = 2, col = "grey50") +
  geom_vline(xintercept = 0, lty = 2, col = "grey50") +
  scale_fill_manual(values = c("gold2", "forestgreen")) +
  scale_shape_manual(values=c(24,21))

# Fig 1C
hull_expr_d <- filter(expr_pca.site_sc, habitat=="dune") %>%
  dplyr::slice(chull(PC1,PC2))
hull_expr_nd <- filter(expr_pca.site_sc, habitat=="non-dune") %>%
  dplyr::slice(chull(PC1,PC2))
plot_1C <-  ggplot(data=expr_pca.site_sc,
                   aes(x = PC1, y = PC2, shape=habitat, fill=habitat)) +
  geom_point(size = 4, alpha=.75) +
  geom_polygon(data=hull_expr_d, alpha=.5) +
  geom_polygon(data=hull_expr_nd, alpha=.5) +
  
  labs(x=paste0("PC1 ","(",round(expr_pve_df[1,2], 2),"%)"),
       y=paste0("PC2 ", "(",round(expr_pve_df[2,2],2),"%)")) +

  theme_bw(base_size = 12) +
  theme(legend.position = "none",
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()) +#,

  geom_hline(yintercept = 0, lty = 2, col = "grey50") +
  geom_vline(xintercept = 0, lty = 2, col = "grey50") +
  scale_fill_manual(values=c("gold2", "forestgreen")) +
  scale_shape_manual(values=c(24,21))

# Fig 1D
hull_rmats_d <- filter(rmats_splice_pca.site_sc, habitat=="dune") %>%
  dplyr::slice(chull(PC1,PC2))
hull_rmats_nd <- filter(rmats_splice_pca.site_sc, habitat=="non-dune") %>%
  dplyr::slice(chull(PC1,PC2))
plot_1D <- ggplot(data=rmats_splice_pca.site_sc,
                                 aes(x = PC1, y = PC2, shape=habitat,
                                     fill=habitat)) +
  geom_point(size = 4, alpha=.75) +
  geom_polygon(data=hull_rmats_d, alpha=0.5) +
  geom_polygon(data=hull_rmats_nd, alpha=0.5) +
  
  labs(x=paste0("PC1 ","(",round(rmats_splice_pve_df[1,2], 2),"%)"),
       y=paste0("PC2 ", "(",round(rmats_splice_pve_df[2,2],2),"%)")) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none",
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()) + 
  
  geom_hline(yintercept = 0, lty = 2, col = "grey50") +
  geom_vline(xintercept = 0, lty = 2, col = "grey50") +
  scale_color_manual(values=c("gold2", "forestgreen")) +
  scale_fill_manual(values=c("gold2", "forestgreen")) +
  scale_shape_manual(values=c(24,21))

#FIG_1_v2 <- gsd_sampling_map + plot_1B + plot_1C + plot_1D + plot_layout(widths = 1)
FIG_1_v2 <- plot_spacer() + plot_1B + plot_1C + plot_1D + plot_layout(widths = 1)
ggsave("figures/FIG_1_v2_raw.pdf", plot=FIG_1_v2,device = "pdf",
       height = 131.25, width = 175, dpi = 300, units = "mm")

#### Figure 2 (Co-expression networks) ####
load("data2/Rdata/WGCNA_networks.mean_counts_3.softPower_18.minModuleSize_30.Rdata.Rdata")
load("data2/Rdata/WGCNA_module_preservation.mean_counts_3.softPower_18.minModuleSize_30.Rdata")

# Panel A
mod_sizes_df <- data.frame(table(d_colors)) %>%
  rename(mod_color=d_colors) %>%
  mutate(Ecotype="Dune") %>%
  full_join(data.frame(table(nd_colors)) %>% 
              rename(mod_color=nd_colors) %>%
              mutate(Ecotype="Non-dune")) %>%
  rename(mod_size=Freq)

jitter <- position_jitter(width = 0.2, height = 0.1)
plot_mod_sizes <- ggplot(data=mod_sizes_df %>% filter(!mod_color %in% c("gold", "grey")),
       aes(x=Ecotype, y=log(mod_size))) +
  geom_point(aes(fill=mod_color), shape=21, size=3, position = jitter) +
  scale_y_continuous(limits = c(3.25,8), breaks=c(4,6,8)) +
  scale_fill_identity() +
  geom_boxplot(alpha=0.5, outlier.shape = NA) +
  theme_bw(base_size = 12) +
  labs(y="log(Module size)")

# Panel B (module preservation)
ref <- 1 #non-dune
test <- 2 #dune 
mod_colors <- rownames(module_preservation$preservation$observed[[ref]][[test]])
mod_sizes <- module_preservation$preservation$Z[[ref]][[test]][, 1];
median_rank <- module_preservation$preservation$observed[[ref]][[test]][, 2]
z_summary <- module_preservation$preservation$Z[[ref]][[test]][, 2]

mod_pres_df <- data.frame(mod_colors, mod_sizes, median_rank, z_summary)

squash_axis <- function(from, to, factor) { 
  # A transformation function that squashes the range of [from, to] by factor on a given axis 
  trans <- function(x) { 
    isq <- x > from & x < to & !is.na(x)
    ito <- x >= to & !is.na(x)
    x[isq] <- from + (x[isq] - from)/factor
    x[ito] <- from + (to - from)/factor + (x[ito] - to)
    return(x)
  }
  inv <- function(x) {
    isq <- x > from & x < from + (to - from)/factor & !is.na(x)
    ito <- x >= from + (to - from)/factor & !is.na(x)
    x[isq] <- from + (x[isq] - from) * factor
    x[ito] <- to + (x[ito] - (from + (to - from)/factor))
    return(x)
  }
  return(trans_new("squash_axis", trans, inv))
}

plot_mod_pres <- ggplot(data=mod_pres_df %>% filter(!mod_colors %in% c("gold", "grey")),
       aes(x=mod_sizes, y=z_summary)) +
  geom_point(aes(fill=mod_colors),shape=21, size=3) +
  scale_fill_identity() +
  theme_bw(base_size = 12) +
  scale_x_continuous(trans = "log10",limits = c(20,2000),
                     breaks = c(20, 100, 500, 2000)) +
  scale_y_continuous(trans = squash_axis(50, 100, 10),
                     limits = c(-2.5,100), breaks = c(10,25,40,50,75,100)) +
  #coord_trans(x = "log10") +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 2, lty=2, col="blue") +
  geom_hline(yintercept = 10, lty=2, col="dark green") +
  labs(x="Module size", y="Preservation Zsummary")
  

plot_mod_sizes + plot_mod_pres
ggsave("figures/FIG_coexpression_raw.pdf", plot=plot_mod_sizes + plot_mod_pres,
       device = "pdf", height = 90, width = 175, dpi = 300, units = "mm")

#### Figure 3 (Manhattan plots) #####
load("data2/Rdata/Fst_2Mb_windows.Rdata") #from analyze_Fst.R
load("data2/Rdata/DESeq2_results_Shrink.Rdata") #from analyze_differential_expression.R
load("data2/Rdata/rmats_results_dfs.Rdata") #from analyze_splicing_rMATS.R
load("data2/Rdata/LFC_dPSI_vs_Fst.Rdata") #from analyze_Fst.R

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

# Panel A
Fst_manhattan_df <- Fst %>% 
  inner_join(cum_lengths, by = "chrom") %>% 
  mutate(cumstart = BIN_START + cumstart) %>%
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

Fst_manhattan_plot <- ggplot(data=Fst_manhattan_df,
                             aes(x=cumstart, y=WEIGHTED_FST)) +
  # alternate shading of chromosomes
  geom_rect(data=cum_lengths, aes(x=NULL,y=NULL,
                                  xmin=cumstart, xmax=cumend,
                                  ymin=0,ymax=.8,
                                  fill=as.factor(chrom)), alpha=0.4) +
  scale_fill_manual(values = rep(c("light grey", "white"), 17 )) +
  
  # add mean Fst line and inversion regions
  geom_line(linewidth=.15) +
  geom_rect(data=inversion_regions, aes(x=NULL, y=NULL,
                                        xmin=inv_cumstart, xmax=inv_cumend,
                                        ymin=0, ymax=0.8),
            fill="red",alpha=0.25) +
  # customize axes
  scale_x_continuous(label = axisdf_Fst$chrom, breaks = axisdf_Fst$center,
                     expand = c(.01,.01)) +
  scale_y_continuous(expand=c(.01,.01)) +
  coord_cartesian(clip = 'off') +
  
  # custom theme
  theme_bw(base_size=12) +
  theme(legend.position = "none",
        #axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),
        #text = element_text(size=24),
        axis.text.x=element_blank(),
        axis.line.x = element_blank(),
        axis.ticks = element_blank(),
        #axis.text.y = element_text(size=18),
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  labs(x="", y=expression(italic(F)["st"]))

# Panel B

transcripts <- data.frame(rtracklayer:: #read-in our annotations, for the gene positions. 
                    import("data/ref_genome_Ha412HO/HAN412_Eugene_curated_v1_1.gff3")) %>%
  rename(chrom=seqnames, Ha412_gene=Parent) %>%
  mutate(Ha412_gene=unlist(Ha412_gene)) %>%
  filter(type %in% c("mRNA", "ncRNA", "rRNA"))

expr_manhattan_df_noLFCthreshold <- left_join(data.frame(de_results_Shrink) %>%
                                                tibble::rownames_to_column("Ha412_gene"),
                                              dplyr::select(transcripts,
                                                            Ha412_gene,
                                                            chrom,
                                                            start,
                                                            end,
                                                            strand)) %>%
  mutate(Ha412_gene=gsub("mRNA","gene",Ha412_gene)) %>%
  mutate(Ha412_gene=gsub("ncRNA","gene",Ha412_gene)) %>%
  inner_join(cum_lengths) %>% #inner join gets rid of non-chromosomal contigs
  mutate(gene_cumstart = start + cumstart, sig=if_else(padj<.05, "yes", "no"))

sig_data_expr <- expr_manhattan_df_noLFCthreshold %>% 
  subset(padj < .05)
notsig_data_expr <- expr_manhattan_df_noLFCthreshold %>% 
  subset(padj >= .05) %>%
  group_by(chrom) %>% 
  sample_frac(0.75)

notsig_remove_expr <- subset(notsig_data_expr,abs(log2FoldChange) < .05) %>%
  sample_frac(0.95) #remove 95% of LFC <.05 non-significant genes.
notsig_keep_expr <- !(notsig_data_expr$Ha412_gene %in% notsig_remove_expr$Ha412_gene)
notsig_data_expr <- notsig_data_expr[notsig_keep_expr,]

expr_manhattan_df_noLFCthreshold_reduced <- bind_rows(sig_data_expr,
                                                      notsig_data_expr)

# get middle position of each chromosome for x axis label
axisdf_expr <- expr_manhattan_df_noLFCthreshold_reduced %>%
  group_by(chrom) %>%
  summarize(center=(max(gene_cumstart) + min(gene_cumstart)) / 2) %>%
  mutate(chrom=c("1","2","3","4","5",
                 "6","7","8","9","10",
                 "11","12","13","14","15",
                 "16","17"))

expr_manhattan_plot <- ggplot(expr_manhattan_df_noLFCthreshold_reduced,
                              aes(x=gene_cumstart, y=log2FoldChange)) +
  #alpha=as.factor(sig)=="yes")) +
  #color=-log10(padj))) +
  # alternate shading of chromosomes
  geom_rect(data=cum_lengths,
            aes(xmin=cumstart, xmax=cumend,
                ymin=min(expr_manhattan_df_noLFCthreshold_reduced$log2FoldChange),
                ymax=max(expr_manhattan_df_noLFCthreshold_reduced$log2FoldChange),
                fill=as.factor(chrom)),
            inherit.aes = F,
            alpha=0.5) +
  scale_fill_manual(values = rep(c("light grey", "white"), 17), guide = "none") +
  
  # mark inversion regions
  geom_rect(data=inversion_regions,
            aes(xmin=inv_cumstart, xmax=inv_cumend,
                ymin=min(expr_manhattan_df_noLFCthreshold_reduced$log2FoldChange),
                ymax=max(expr_manhattan_df_noLFCthreshold_reduced$log2FoldChange)),
            fill="red", alpha=0.25, inherit.aes = F) +
  
  # add each gene
  geom_point(aes(size=as.factor(sig=="yes")), shape=1, stroke = .25) + #alpha=0.75, aes(size=-log10(padj))) +
  #scale_alpha_manual(name="FDR < .05",values = c(.1,1), labels=c("no", "yes")) +
  scale_size_manual(name="FDR < .05", values = c(.5,2), labels=c("no", "yes")) +
  #scale_color_gradientn(colors=c("grey","turquoise","blue","navyblue","black")) +
  
  # custom axes
  scale_x_continuous(label = axisdf_expr$chrom,
                     breaks = axisdf_expr$center, expand=c(.01,.01)) +
  scale_y_continuous(expand=c(.01,.01)) +
  coord_cartesian(clip = 'off') +
  
  # custom theme
  theme_bw(base_size=12) +
  theme(legend.position = "none",
        #axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),
        #text = element_text(size=24),
        axis.line.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x=element_blank(),
        #plot.title = element_text(size=18),
        #axis.text.y = element_text(size=18),
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  
  labs(x="", y="log2FC")

# Panel C
splice_manhattan_df <- inner_join(all_AS_events_deltaPSI, cum_lengths) %>%
  #rowwise() %>%
  mutate(event_startcum = event_start + cumstart, sig=if_else(FDR<.05, "yes", "no"))
axisdf_splice <- splice_manhattan_df %>% #midpoint for each chrom
  group_by(chrom) %>%
  summarize(center=(max(event_startcum) + min(event_startcum)) / 2) %>%
  mutate(chrom=c("1","2","3","4","5",
                 "6","7","8","9","10",
                 "11","12","13","14","15",
                 "16","17"))
# split apart significant DE genes and nonsignifcant DE genes 
# so we can downsample the nonsig data, for better vis
sig_splice_data <- splice_manhattan_df %>% 
  subset(FDR < 0.05)
notsig_splice_data <- splice_manhattan_df %>% 
  subset(FDR >= 0.05) %>%
  group_by(chrom) %>% 
  sample_frac(0.5)
splice_manhattan_df_reduced <- bind_rows(sig_splice_data, notsig_splice_data) 

splice_manhattan_plot <- ggplot(splice_manhattan_df_reduced,
                                aes(x=event_startcum, y=IncLevelDifference),
                                alpha=IncLevelDifference) +
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
  geom_point(aes(size=as.factor(sig=="no")), shape=1, stroke = .25) +
  #scale_alpha_manual(name="FDR < .05",values = c(.1,.5), labels=c("no", "yes")) +
  scale_size_manual(name="Significance (FDR < .05)", values = c(2,.5), labels=c("yes", "no")) +
  #geom_point(alpha=0.75, aes(size=-log(FDR))) +
  #scale_color_gradientn(colors=c("grey","turquoise","blue","navyblue","black")) +
  
  # custom axes
  scale_x_continuous(label = axisdf_splice$chrom,
                     breaks = axisdf_splice$center, expand=c(0.01,0.01)) +
  scale_y_continuous(expand=c(0.01,0.01)) +
  
  # custom theme
  theme_bw(base_size = 12) +
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
        legend.text = element_text(size=10),
        legend.title = element_text(size=10),
        #axis.text = element_text(size=18),
        legend.box.margin = margin(t = -15)) +
  
  labs(x="Chromosome",y=~paste(Delta,"PSI"))

# Panel D
plot_LFC_vs_Fst <- ggplot(data=LFC_vs_Fst_df, aes(x=Fst, y=abs(log2FoldChange))) +
  geom_point(size=1.5, alpha=.25, stroke=.125) +
  geom_smooth(method = "lm", color="red", linewidth=.75) +
  theme_bw(base_size = 12) +
  #ylim(c(0,10)) +
  theme(#text = element_text(size=24),#,
    #legend.title=element_blank(),
    #plot.title = element_text(size=18),
    #axis.text = element_text(size=18),
    panel.grid.major = element_blank()) +
  labs(x=expression(italic(F)["st"]), y="|log2FC|")

# Panel E
plot_dPSI_vs_Fst <- ggplot(data=dPSI_vs_Fst_df, aes(x=Fst, y=max_deltaPSI)) +
  geom_point(size=1.5, alpha=.25, stroke=.125) +
  geom_smooth(method = "lm", color="red", linewidth=.75) +
  theme_bw(base_size=12) +
  theme(#text = element_text(size=24),#,
    #legend.title=element_blank(),
    #plot.title = element_text(size=18),
    #axis.text = element_text(size=18),
    panel.grid.major = element_blank()) +
  labs(x=expression(italic(F)["st"]), y=~paste("|",Delta,"PSI|"))

plot_manhattans <- (Fst_manhattan_plot / expr_manhattan_plot / splice_manhattan_plot) /
  (plot_spacer() | plot_LFC_vs_Fst | plot_dPSI_vs_Fst) +
  plot_layout(heights = c(1,1,1,1.25))

ggsave("figures/plot_manhattans_raw.pdf", plot = plot_manhattans,
       device = "pdf", width = 175, height = 175, units = "mm", dpi = 300)

#### Figure 4 (GO enrichment) ####
load("data2/Rdata/GO_results_main_text.Rdata")

# For this figure we will save each panel individually
# and then combine with Inkscape.

# Panel A
GO_DE_dune_dotplot <- ggplot(data = GO_results_DE_dune_ALL,
                             aes(x=fold_enrichment,
                                 y=fct_reorder(Description, -log10(p.value),
                                               .desc=F))) +
  geom_point(aes(size=-log10(p.value), color=Type)) +
  scale_color_manual(name = "GO category", values = c("black", "grey50", "grey75"),
                     labels = c("Biological process",
                                "Cellular component",
                                "Molecular function"), guide = "none") +
  scale_size_continuous(name = bquote(paste("-",log[10], "(p-value)")),
                        range = c(1,4), breaks = c(2,4,8,16) ) + 
  
  scale_x_continuous(expand = c(.2,0), breaks = c(2,4,6,8)) +
  scale_y_discrete(expand = c(0.1,.1)) +
  coord_fixed(1.5) +
  theme_bw(base_size = 12) +
  theme(axis.text.y = element_text(size=8),
        legend.text = element_text(size=8),
        legend.position = "none",
        legend.direction = "vertical",
        legend.title = element_text(size=8),
        legend.margin=margin(t=-10)) +

  guides(color = guide_legend(override.aes = list(size=3))) +
  labs(y="", x="")
ggsave("figures/plot_GO_DE_dune_raw.pdf", plot = GO_DE_dune_dotplot,
       device = "pdf", width = 87.5, dpi = 300, units = "mm")


# Panel B
venn_dia_DE_DS <- euler(c(DE=4871, DS=806, "DE&DS"=232))
d <- eulerr:::plot.euler(venn_dia_DE_DS,
                         labels = NA,
                         quantities = NA,
                         fills = list(fill = c("#317EC2", "red", "purple"),
                                      alpha=0.5))
pdf(file = "figures/venn_dia_DE_DS_blank.pdf", width = 6, height = 4.5)#, bg = "transparent")
d
dev.off()


# Panel C
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
                        range = c(1,6), breaks = c(2,4,8,16)) + 
  
  scale_x_continuous(expand = c(.2,0)) +
  scale_y_discrete(expand = c(.075,0)) +
  coord_fixed(4) +
  
  theme_bw(base_size=12) +
  theme(axis.title.x = element_text(margin=margin(l=0,r=200)),
        axis.text.y = element_text(size=8),
        legend.text = element_text(size=8),
        legend.position = "bottom",
        legend.direction = "vertical",
        legend.title = element_text(size=8),
        legend.margin=margin(t=-10)) +
  guides(color = guide_legend(override.aes = list(size=3))) +
  labs(x="Fold enrichment", y="")
ggsave("figures/plot_GO_DS_raw.pdf", plot = GO_DS_dotplot,
       device = "pdf", width = 87.5, dpi = 300, units = "mm")


# Panel D
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
                        range = c(1,4.5), breaks = c(2,4,8,16)) + 
  
  scale_x_continuous(expand = c(.1,0)) +
  scale_y_discrete(expand = c(.25,0), position = "right") +
  coord_fixed(1) +
  
  theme_bw(base_size = 12) +
  theme(axis.text.y = element_text(size=8),
        legend.text = element_text(size=8),
        legend.position = "none",
        legend.direction = "vertical",
        legend.title = element_text(size=8),
        legend.margin=margin(t=-10)) +  
  guides(color = guide_legend(override.aes = list(size=3))) +
  labs(x="", y="")
ggsave("figures/plot_GO_intersect_raw.pdf", plot = GO_intersect_dotplot,
       device = "pdf", width = 87.5, dpi = 300, units = "mm")

#### Figure 5. GLH17 ####
load("data2/Rdata/rmats_results_dfs.Rdata") #from analyze_splicing_rMATS.R
load("data2/Rdata/DESeq2_normalized_counts.Rdata") #from analyze differential splicing

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
  mutate(ecotype=if_else(sample %in% dune_fq, "Dune", "Non-dune"),
         PSI=as.numeric(PSI))

norm_expr_df <- normalized_counts[,order(colnames(normalized_counts))] %>% 
  rownames_to_column("Ha412_gene") %>%
  mutate(Ha412_gene=gsub(".*RNA","gene",Ha412_gene)) %>%
  pivot_longer(2:25, names_to = "sample", values_to = "normalized_expression") %>%
  mutate(ecotype=if_else(grepl("non", sample), "Non-dune", "Dune"),
         Ha412_gene=gsub("gene:","",Ha412_gene))

plot_GLH17_psi <- ggplot(data=PSI_df %>% filter(ID=="SE_8901"),
                         aes(x=ecotype, y=PSI,
                             fill=ecotype, shape=ecotype)) +
  geom_point(size=3, alpha=.75, position = position_jitter(width = .2)) +
  geom_boxplot(alpha=0.25, outlier.shape = NA, linewidth=.25) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("gold2", "forestgreen")) +
  scale_color_manual(values = c("gold2", "forestgreen")) +
  scale_shape_manual(values=c(24,21)) +
  labs(x="Ecotype", y="Percent Spliced In")

plot_GLH17_expr <- ggplot(data=norm_expr_df %>%
                            filter(Ha412_gene %in% c("Ha412HOChr02g0088581")),
                          aes(x=ecotype, y=normalized_expression,
                              shape=ecotype, fill=ecotype)) +
  geom_point(size=3, alpha=.75, position = position_jitter(width = .2)) +
  geom_boxplot(alpha=0.25, outlier.shape = NA, linewidth=.25) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("gold2", "forestgreen")) +
  scale_color_manual(values = c("gold2", "forestgreen")) +
  scale_shape_manual(values=c(24,21)) +
  labs(x="Ecotype", y="Norm. expression")

plot_GLH17_expr + plot_GLH17_psi
ggsave(filename = "figures/plot_GLH17_raw.pdf", plot = plot_GLH17_expr + plot_GLH17_psi,
       device = "pdf", dpi = 300, width = 87.5, height = 65.625, units = "mm")
