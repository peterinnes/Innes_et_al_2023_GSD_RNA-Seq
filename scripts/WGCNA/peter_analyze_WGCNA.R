# This script performs gene coexpression network analysis for dune versus
# non-dune samples.
# by Kaylee Rosenberger and Peter Innes

library(DESeq2)
library(WGCNA)
library(dplyr)

#### 1. data processing and prep ####
# read-in, filter, and transform raw count data
sample_table_deseq2 <- read.table("analysis/DESeq2/sample_table.txt", header=T)
sample_table_deseq2$habitat <- as.factor(sample_table_deseq2$habitat)
counts <- read.table("analysis/DESeq2/htseq-count_out/htseq-count_results.2022-6-26.txt",
                     row.names = 1)

names(counts) <- paste(sample_table_deseq2$habitat,"_",
                       sample_table_deseq2$sample_no, sep = "")

# filter the raw counts, mean counts per sample
keep_wgcna <- which(rowMeans(as.matrix(counts)) >= 1) #same prefilter as for DESeq2
keep_stringent_wgcna <- which(rowMeans(as.matrix(counts)) >= 3) #more stringent
keep_vstringent_wgcna <- which(rowMeans(as.matrix(counts)) >= 10) #very stringent, for iterWGCNA

filtered_counts_wgcna <- counts[keep_wgcna,] %>%
  head(-3) #remove the last three rows, "__no_feature", "__ambiguous" "__alignment_not_unique"
filtered_stringent_counts_wgcna <- counts[keep_stringent_wgcna,] %>%
  head(-3)
filtered_vstringent_counts_wgcna <- counts[keep_vstringent_wgcna,] %>%
  head(-3) %>% 
  mutate(d_n_zero=rowSums(dplyr::select(.,starts_with("dune")) == 0),
         nd_n_zero=rowSums(dplyr::select(.,starts_with("non-dune")) == 0)) %>%
  filter(d_n_zero <= 6 & nd_n_zero <= 6) %>%
  dplyr::select(!c(d_n_zero, nd_n_zero))
  
dim(filtered_counts_wgcna) #32308
dim(filtered_stringent_counts_wgcna) #28770
dim(filtered_vstringent_counts_wgcna) #24421

rlog_filtered_counts_wgcna <- rlog(as.matrix(filtered_stringent_counts_wgcna)) %>%
  data.frame()
rlog_filtered_vstringent_counts_wgcna <- rlog(as.matrix(filtered_vstringent_counts_wgcna)) %>%
  data.frame()

write.table(rlog_filtered_counts_wgcna,
            file="analysis/WGCNA/rlog_filtered_counts_wgcna.tsv",
            quote = F, row.names = T, sep = "\t", col.names = T)
write.table(rlog_filtered_vstringent_counts_wgcna,
            file="analysis/WGCNA/rlog_filtered_vstringent_counts_wgcna.tsv",
            quote = F, row.names = T, sep = "\t", col.names = T)

#### write vstringent expr data for iterativeWGCNA ####
nd_expr_vstringent <- data.frame(rlog_filtered_vstringent_counts_wgcna) %>%
  dplyr::select(starts_with("non.dune")) %>%
  rownames_to_column("gene")
write.table(x = nd_expr_vstringent,
            file = "analysis/WGCNA/iterativeWGCNA/non_dune_expr_data.iterWGCNA.tsv",
            quote = F, sep="\t", row.names = F)
d_expr_vstringent <- data.frame(rlog_filtered_vstringent_counts_wgcna) %>%
  dplyr::select(starts_with("dune")) %>%
  rownames_to_column("gene")
write.table(x = d_expr_vstringent,
            file = "analysis/WGCNA/iterativeWGCNA/dune_expr_data.iterWGCNA.tsv",
            quote = F, sep="\t", row.names = F)

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

# check scale-free network topology
k1=softConnectivity(datExpr = multi_expr[[1]]$data, power=18) #non-dune
k2=softConnectivity(datExpr = multi_expr[[2]]$data, power=18) #dune
pdf(file="figures/scale_free_topology.pdf", width=3.45, height=1.75)
par(mfrow=c(1,2))
scaleFreePlot(k1)
scaleFreePlot(k2)

#### 3. build networks and plot dendrograms ####
nd_network <- blockwiseModules(multi_expr[[1]]$data, maxBlockSize = 30000,
                              power = 18, minModuleSize = 30,
                              networkType = "signed", deepSplit = 2,
                              TOMType = "signed",mergeCutHeight = 0.25,
                              numericLabels = TRUE, minKMEtoStay = 0.3,
                              saveTOMs = FALSE,  verbose = 4, nThreads = 20)

d_network <- blockwiseModules(multi_expr[[2]]$data, maxBlockSize = 30000,
                              power = 18, minModuleSize = 30,
                              networkType = "signed", deepSplit = 2,
                              TOMType = "signed",mergeCutHeight = 0.25,
                              numericLabels = TRUE, minKMEtoStay = 0.3,
                              saveTOMs = FALSE, verbose = 4, nThreads = 20)

# module labels as colors for nondune set
nd_colors = labels2colors(nd_network$colors)
# show number and size of modules
table(nd_network$colors)

# module labels as colors for dune set
d_colors = labels2colors(d_network$colors)
# matching labels to those used in non-dune set 
d_colors = matchLabels(d_colors, nd_colors)
# show number and size of modules
table(d_network$colors)

# make a dataframe of the module colors and sizes
mod_sizes_df <- data.frame(table(d_colors)) %>%
  dplyr::rename(mod_color=d_colors) %>%
  mutate(Ecotype="Dune") %>%
  full_join(data.frame(table(nd_colors)) %>% 
              dplyr::rename(mod_color=nd_colors) %>%
              mutate(Ecotype="Non-dune")) %>%
  dplyr::rename(mod_size=Freq)

# exclude gold (module created with random genes) and grey (module of all unplaced genes) 
mod_sizes_df %>% filter(!mod_color %in% c("gold", "grey")) %>%
  group_by(Ecotype) %>%
  summarise(Mean=mean(mod_size), N=n())

fit_mod_sizes <- lm(log(mod_size) ~ Ecotype, data=mod_sizes_df)
summary(fit_mod_sizes)

## average module size when not excluding grey and gold
#mean(table(nd_network$colors))
#mean(table(d_network$colors))

# save results 
save(nd_network, nd_colors, d_network, d_colors,
     file = "data2/Rdata/WGCNA_networks.mean_counts_3.softPower_18.minModuleSize_30.Rdata")

#save(nd_network, file = "analysis/WGCNA/non-dune_wgcna_network.Rdata")
#save(d_network, file= "analysis/WGCNA/dune_wgcna_network.Rdata")

# plot non-dune dendrogram
png(filename = "figures/plot_nondune_dendrogram.png", width = 250,
    height = 75, res = 300, units = "mm")
plotDendroAndColors(nd_network$dendrograms[[1]], cbind(nd_colors, d_colors),
                    c("Non-dune", "Dune"),
                    main = "Non-dune network",
                    dendroLabels = FALSE, cex.colorLabels = 1,
                    cex.axis=1, cex.lab=1, cex.main=1)
dev.off()

# plot dune dendrogram with module labels matched to non-dune labels
png(filename = "figures/plot_dune_dendrogram.png", width = 250,
    height = 75, res = 300, units = "mm")
plotDendroAndColors(d_network$dendrograms[[1]], cbind(d_colors, nd_colors),
                    c("Dune", "Non-dune"),
                    main = "Dune network",
                    dendroLabels = FALSE, cex.colorLabels = 1,
                    cex.axis=1, cex.lab=1, cex.main=1)

dev.off()

#### 4. analyze module preservation ####

# number of modules with the same labels 
length(intersect(unique(d_colors), unique(nd_colors)))

# function for unique modules
outersect <- function(x, y) {
  sort(c(setdiff(x, y),
         setdiff(y, x)))
}

outersect(unique(d_colors), unique(nd_colors))

# Determining unique modules to each network: 
# modules unique to non-dune
nd_unique_mods = setdiff(unique(nd_colors), unique(d_colors))

# modules unique to dune 
d_unique_mods = setdiff(unique(d_colors), unique(nd_colors))

# calculate module preservation
multi_color = list("Non-dune" = nd_colors, "Dune" = d_colors)
module_preservation = modulePreservation(multi_expr, multi_color,
                                         dataIsExpr = TRUE, 
                                         networkType = "signed",
                                         referenceNetworks = c(1:2),
                                         nPermutations = 200, randomSeed = 1,
                                         quickCor = 0, verbose = 3,
                                         maxModuleSize = 2000)

# save results because it takes a long time to run 
save(module_preservation, file = "data2/Rdata/WGCNA_module_preservation.mean_counts_3.softPower_18.minModuleSize_30.Rdata")

# plotting mod preservation ####
# How well preserved are non-dune modules in the dune network?

ref <- 1 #non-dune
test <- 2 #dune 
mod_colors = rownames(module_preservation$preservation$observed[[ref]][[test]])
mod_sizes = module_preservation$preservation$Z[[ref]][[test]][, 1];

# leave out grey (unplaced genes) and gold (module built from random genes)
plot_mods = !(mod_colors %in% c("grey", "gold"));

# Text labels for points
text = mod_colors[plot_mods];

# Auxiliary convenience variable
plot_data_mod_pres = cbind(module_preservation$preservation$observed[[ref]][[test]][, 2], module_preservation$preservation$Z[[ref]][[test]][, 2])

# Main titles for the plot
mains = c("Preservation Median rank", "Preservation Zsummary");

# Start the plot
png(filename="figures/module_preservation.png", width=12, height=4.5, res=300,
    units="in")

par(mfrow = c(1,2))
par(mar = c(4.5,4.5,2.5,2.5))
for(p in 1:2) {
  min = min(plot_data_mod_pres[, p], na.rm = TRUE);
  max = max(plot_data_mod_pres[, p], na.rm = TRUE);
  # Adjust ploting ranges appropriately
  if (p==2)
  {
    if (min > -max/10) min = -max/10
    ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
  } else
    ylim = c(max + 0.1 * (max-min), min - 0.1 * (max-min))
  plot(mod_sizes[plot_mods], plot_data_mod_pres[plot_mods, p], col = 1,
       bg = mod_colors[plot_mods], pch = 21,
       cex=2.25,
       ylab = mains[p], xlab = "Module size", log = "x",
       ylim = ylim,
       xlim = c(20, 2000), cex.lab = 1.75, cex.axis = 1.25, cex.main = 1.75)
  #labelPoints(mod_sizes[plotMods], plotData[plotMods, p], text, cex = 1, offs = 0.08)
  # For Zsummary, add threshold lines
  if (p==2)
  {
    abline(h=0, lwd=2)
    abline(h=2, col = "blue", lty = 2, lwd=2)
    abline(h=10, col = "darkgreen", lty = 2, lwd=2)
  }
}
dev.off()



