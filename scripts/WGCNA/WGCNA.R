# This script builds gene coexpression networks for dune versus
# non-dune samples and tests preservation of non-dune modules in the dune network
# by Kaylee Rosenberger and Peter Innes
# Most of the code adapted from WGCNA tutorials (Langfelder et al)

library(DESeq2)
library(WGCNA)
library(dplyr)
library(ggplot2)
library(scales)
library(patchwork)

#### 1. data processing and prep ####
# read-in, filter, and transform raw count data
sample_table_deseq2 <- read.table("data/sample_table.txt", header=T)
sample_table_deseq2$habitat <- as.factor(sample_table_deseq2$habitat)
counts <- read.table("data/htseq-count_results.2022-6-26.txt",
                     row.names = 1)

names(counts) <- paste(sample_table_deseq2$habitat,"_",
                       sample_table_deseq2$sample_no, sep = "")

# filter the raw counts, mean counts per sample. very stringent filter
keep_vstringent_wgcna <- which(rowMeans(as.matrix(counts)) >= 10) 

filtered_vstringent_counts_wgcna <- counts[keep_vstringent_wgcna,] %>%
  head(-3) %>% 
  mutate(d_n_zero=rowSums(dplyr::select(.,starts_with("dune")) == 0),
         nd_n_zero=rowSums(dplyr::select(.,starts_with("non-dune")) == 0)) %>%
  filter(d_n_zero <= 6 & nd_n_zero <= 6) %>%
  dplyr::select(!c(d_n_zero, nd_n_zero))
  
dim(filtered_vstringent_counts_wgcna) #24421

rlog_filtered_vstringent_counts_wgcna <- DESeq2::rlog(as.matrix(filtered_vstringent_counts_wgcna)) %>%
  data.frame()

#### 2. reformat expression data into a list of lists for WGCNA check data ####
non_dune_expr_data <- data.frame(t(rlog_filtered_vstringent_counts_wgcna %>%
  dplyr::select(starts_with("non.dune"))))
dune_expr_data <- data.frame(t(rlog_filtered_vstringent_counts_wgcna %>% 
  dplyr::select(starts_with("dune"))))

num_sets <- 2 #two sets, dune and non dune
set_labels <- c("Non-dune", "Dune")
multi_expr <- list(non_dune = list(data = non_dune_expr_data),
                   dune = list(data = dune_expr_data))
                   
checkSets(multi_expr) #check that formatting is correct
# check for genes with too many missing samps. All look OK.
goodSamplesGenesMS(multi_expr, verbose = 3)

# Check clustering of samples, based on euclidian distance,
# to determine if there are any outliers
sample_trees = list()
for (set in 1:num_sets) {
  sample_trees[[set]] = hclust(dist(multi_expr[[set]]$data), method = "average")
}

# plotting cluster dendrograms of samples 
par(mfrow=c(2,1))
par(mar = c(0, 4, 2, 0))
for(set in 1:num_sets) {
  plot(sample_trees[[set]], main = paste("Sample clustering on all genes in", set_labels[set]),
       xlab="", sub="", cex = 0.7)
}

#### 2.5 scale-free topology and determine soft thresholding power ####
# NOTE: Only need to run this chunk of code one time
# to determine an appropriate soft-thresholding power 

# choose a set of soft-thresholding powers to try (default from tutorial)
powers = c(seq(4,10,by=1), seq(12,30, by=2))
power_tables = vector(mode = "list", length = num_sets)
for(set in 1:num_sets) {
  power_tables[[set]] = list(data = pickSoftThreshold(multi_expr[[set]]$data, powerVector=powers, 
                                                     verbose = 2, networkType = "signed")[[2]])
}
collectGarbage()

# Plotting results 
colors = c("black", "red")

plot_cols = c(2,5,6,7) #plot these columns of the returned scale free analysis tables
col_names = c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity",
             "Max connectivity");

# Get the minima and maxima of the plotted points
ylim = matrix(NA, nrow = 2, ncol = 4);
for (set in 1:num_sets)
{
  for (col in 1:length(plot_cols))
  {
    ylim[1, col] = min(ylim[1, col], power_tables[[set]]$data[, plot_cols[col]], na.rm = TRUE);
    ylim[2, col] = max(ylim[2, col], power_tables[[set]]$data[, plot_cols[col]], na.rm = TRUE);
  }
}
# Plot the quantities in the chosen columns vs. the soft thresholding power
par(mfcol = c(2,2));
par(mar = c(4.2, 4.2 , 2.2, 0.5))
cex1 = 0.7;
for (col in 1:length(plot_cols)) for (set in 1:num_sets)
{
  if (set==1)
  {
    plot(power_tables[[set]]$data[,1], -sign(power_tables[[set]]$data[,3])*power_tables[[set]]$data[,2],
         xlab="Soft Threshold (power)", ylab=col_names[col],
         type="n", ylim = ylim[, col],
         main = col_names[col]) + abline(h=0.8, col="red")
    addGrid();
  }
  if (col==1)
  {
    text(power_tables[[set]]$data[,1], -sign(power_tables[[set]]$data[,3])*power_tables[[set]]$data[,2],
         labels=powers,cex=cex1,col=colors[set])
  } else
    text(power_tables[[set]]$data[,1], power_tables[[set]]$data[,plot_cols[col]],
         labels=powers,cex=cex1,col=colors[set])
  if (col==1)
  {
    legend("bottomright", legend = set_labels, col = colors, pch = 20) 
  } else
    legend("topright", legend = set_labels, col = colors, pch = 20) 
}

# Choose power 18 for both based on plots above Choose first power that 
# achieves a good model fit of scale-free topology (first power with an R^2 ~ 0.8)
# I chose the same value for both networks to more easily compare them, and because 18 is a
# pretty good value for both. Also, 18 is recommended by the WGCNA developers for
# smaller sample sizes: 
# https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/faq.html

# check scale-free network topology. for sup mat
k1=softConnectivity(datExpr = multi_expr[[1]]$data, power=18) #non-dune
k2=softConnectivity(datExpr = multi_expr[[2]]$data, power=18) #dune
scaleFreePlot(k1)
scaleFreePlot(k2)

#### 3. build networks and plot dendrograms ####
nd_network <- blockwiseModules(multi_expr[[1]]$data, maxBlockSize = 25000,
                              power = 18, minModuleSize = 30, 
                              mergeCutHeight = 0.25,
                              networkType = "signed", TOMType = "signed",
                              numericLabels = TRUE, saveTOMs = FALSE,
                              verbose = 4, nThreads = 8)

d_network <- blockwiseModules(multi_expr[[2]]$data, maxBlockSize = 25000,
                              power = 18, minModuleSize = 30,
                              mergeCutHeight = 0.25,
                              networkType = "signed", TOMType = "signed",
                              numericLabels = TRUE, saveTOMs = FALSE,
                              verbose = 4, nThreads = 8)

# module labels as colors for nondune set
nd_colors = labels2colors(nd_network$colors)
# show number and size of modules
table(nd_network$colors)

# module labels as colors for dune set
d_colors = labels2colors(d_network$colors)
# matching labels to those used in non-dune set. 
# This step not actually required for any of our reported results. 
d_colors = matchLabels(d_colors, nd_colors)

# make a dataframe of the module colors and sizes
mod_sizes_df <- data.frame(table(d_colors)) %>%
  dplyr::rename(mod_color=d_colors) %>%
  mutate(Ecotype="Dune") %>%
  full_join(data.frame(table(nd_colors)) %>% 
              dplyr::rename(mod_color=nd_colors) %>%
              mutate(Ecotype="Non-dune")) %>%
  dplyr::rename(mod_size=Freq)

# exclude gold (module created with random genes) 
# and grey (module of all unplaced genes) 
mod_sizes_df %>% filter(!mod_color %in% c("gold", "grey")) %>%
  group_by(Ecotype) %>%
  summarise(Mean=mean(mod_size), N=n())

fit_mod_sizes <- lm(mod_size ~ Ecotype,data=mod_sizes_df %>%
                        filter(!mod_color %in% c("gold", "grey")))
summary(fit_mod_sizes) # no significant difference

# plot non-dune dendrogram for supp mat
png(filename = "figures/nondune_dendrogram_R1.png", width = 250,
    height = 75, res = 300, units = "mm")
plotDendroAndColors(nd_network$dendrograms[[1]], nd_colors,
                    "Non-dune", main = "Non-dune dendrogram and module labels",
                    dendroLabels = FALSE, cex.colorLabels = 1,
                    cex.axis=1, cex.lab=1, cex.main=1)
dev.off()

# plot dune dendrogram with non-dune modules labels for supp mat
png(filename = "figures/dune_dendrogram_R1.png", width = 250,
    height = 75, res = 300, units = "mm")
plotDendroAndColors(d_network$dendrograms[[1]], nd_colors,
                    "Non-dune", 
                    main = "Dune dendrogram with non-dune module labels",
                    dendroLabels = FALSE, cex.colorLabels = 1,
                    cex.axis=1, cex.lab=1, cex.main=1)

dev.off()

#### 4. analyze module preservation ####

## number of modules with the same labels 
#length(intersect(unique(d_colors), unique(nd_colors)))
#
## function for unique modules
#outersect <- function(x, y) {
#  sort(c(setdiff(x, y),
#         setdiff(y, x)))
#}
#
#outersect(unique(d_colors), unique(nd_colors))
#
## Determining unique modules to each network: 
## modules unique to non-dune
#nd_unique_mods = setdiff(unique(nd_colors), unique(d_colors))
#
## modules unique to dune 
#d_unique_mods = setdiff(unique(d_colors), unique(nd_colors))

# calculate module preservation for non-dune modules in dune data
multi_color = list(non_dune = nd_colors)
module_preservation = modulePreservation(multi_expr, multi_color,
                                         dataIsExpr = TRUE, 
                                         networkType = "signed",
                                         referenceNetworks = 1,
                                         nPermutations = 100, randomSeed = 1,
                                         quickCor = 0, verbose = 3,
                                         maxModuleSize = 2000)

save(module_preservation, file = "data/vanillaWGCNA_module_preservation.Rdata")
#load("data/iterWGCNA_module_preservation.Rdata")
#load("data/vanillaWGCNA_module_preservation.Rdata")

#### plotting mod preservation ####
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

plot_mod_pres <- ggplot(data=mod_pres_df %>%
                            filter(!mod_colors %in% c("gold", "grey")),
                        aes(x=mod_sizes, y=z_summary)) +
    geom_point(aes(fill=mod_colors),shape=21, size=2) +
    scale_fill_identity() +
    theme_bw(base_size = 8) +
    scale_x_continuous(trans = "log10",limits = c(20,2000),
                       breaks = c(20, 100, 500, 2000)) +
    scale_y_continuous(trans = squash_axis(50, 100, 10),
                       limits = c(-2.5,100), breaks = c(10,25,40,50,75,100)) +
    #coord_trans(x = "log10") +
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = 2, lty=2, col="blue") +
    geom_hline(yintercept = 10, lty=2, col="dark green") +
    labs(x="Module size", y="Preservation Zsummary")

plot_mod_pres_violin <- ggplot(data=mod_pres_df %>%
                            filter(!mod_colors %in% c("gold", "grey")) %>% 
                            mutate(Comparison=""),
                            aes(x=Comparison, y=z_summary)) + 
    geom_violin(alpha=0.25, fill="grey") +
    geom_boxplot(width=.2) +
    geom_hline(yintercept = 10, lty=2, color="dark green") +
    geom_hline(yintercept = 2, lty=2, color="blue") +
    scale_y_continuous(limits = c(-1.25,100)) +
    theme_bw(base_size = 8) +
    labs(x = "Non-dune vs Dune", y = "Preservation Zsummary")

ggsave("figures/FIG_coexpression_R1.pdf", plot = plot_mod_pres_violin + plot_mod_pres,
       device = "pdf", height = 57.75, width = 115.5, dpi = 300, units = "mm")
