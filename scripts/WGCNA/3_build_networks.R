################################################################################################
#WGCNA (weighted gene co-expression network analysis) for Kane lab rotation project
#By Kaylee Rosenberger, following WGCNA package tutorial 
#using Great Sand Dunes National Park prairie sunflower genome (from Peter Innes)
#dataset is a table of 30000+ genes for 25 samples of dune and non-dune populations 
#each cell represents regularized log transformed expression counts 

#This script runs module preservation analyses to determine differences between the modules 
#found in the dune and nondune datasets 

#import libraries
library(BiocManager)
library(WGCNA)
library(tibble)
library(dplyr)

##########################################################################################
#load in data
load("data/Rdata/multi_expr.Rdata")
#important option
options(stringsAsFactor = FALSE)
disableWGCNAThreads()

##########################################################################################
#Scale-free toplogy analysis and choosing soft-thresholding power 
#NOTE: Only need to run this chunk of code one time to determine an appropriate soft-thresholding power 

#choose a set of soft-thresholding powers to try (default from tutorial)
powers = c(seq(4,10,by=1), seq(12,30, by=2))
powerTables = vector(mode = "list", length = num_sets)
for(set in 1:num_sets) {
  powerTables[[set]] = list(data = pickSoftThreshold(multi_exp[[set]]$data, powerVector=powers, 
                                                     verbose = 2, networkType = "signed")[[2]])
}
collectGarbage()

#Plotting results 
colors = c("black", "red")

# Will plot these columns of the returned scale free analysis tables
plotCols = c(2,5,6,7)
colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity",
             "Max connectivity");

# Get the minima and maxima of the plotted points
ylim = matrix(NA, nrow = 2, ncol = 4);
for (set in 1:num_sets)
{
  for (col in 1:length(plotCols))
  {
    ylim[1, col] = min(ylim[1, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
    ylim[2, col] = max(ylim[2, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
  }
}
# Plot the quantities in the chosen columns vs. the soft thresholding power
pdf(file = "plots/assumptions/scaleFreeAnalysis.pdf", wi = 8, he = 6);
par(mfcol = c(2,2));
par(mar = c(4.2, 4.2 , 2.2, 0.5))
cex1 = 0.7;
for (col in 1:length(plotCols)) for (set in 1:num_sets)
{
  if (set==1)
  {
    plot(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
         xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = ylim[, col],
         main = colNames[col]) + abline(h=0.8, col="red")
    addGrid();
  }
  if (col==1)
  {
    text(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
         labels=powers,cex=cex1,col=colors[set])
  } else
    text(powerTables[[set]]$data[,1], powerTables[[set]]$data[,plotCols[col]],
         labels=powers,cex=cex1,col=colors[set])
  if (col==1)
  {
    legend("bottomright", legend = set_labels, col = colors, pch = 20) 
  } else
    legend("topright", legend = set_labels, col = colors, pch = 20) 
}
dev.off()

#Choose power 16 for both based on plot scaleFreeAnalysis_mp.pdf--Choose first power that 
#achieves a good model fit of scale-free topology (first power with an R^2 above 0.8)
#I chose the same value for both networks to more easily compare them, and because 16 is a
#pretty good value for both 

##########################################################################################
#Create signed, weighted gene co-expression networks for dune and nondune independently to 
#determine how modules are conserved
#Here we use a soft-thresholding power of 16. Use a max block size of 35000 (> # genes in dataset) 
#so that modules are created in a single block. Use a larger minimum module size of 100
#to more easily compare modules between sets, and a moderate sensitivity of splitting modules with 
#deepSplit = 2. Modules with a height of < 0.25 are merged into one 

#nondune construction
nondune_network = blockwiseModules(multi_exp[[1]]$data, maxBlockSize = 35000, power = 16, 
                                    minModuleSize = 100, networkType = "signed",
                                    deepSplit = 2, TOMType = "signed",
                                    mergeCutHeight = 0.25, numericLabels = TRUE,
                                    minKMEtoStay = 0,
                                    saveTOMs = FALSE,
                                    verbose = 4)

#dune construction
dune_network = blockwiseModules(multi_exp[[2]]$data, maxBlockSize = 35000, power = 16,
                                   minModuleSize = 100, networkType = "signed",
                                   deepSplit = 2, TOMType = "signed",
                                   mergeCutHeight = 0.25, numericLabels = TRUE,
                                   minKMEtoStay = 0,
                                   saveTOMs = FALSE, 
                                   verbose = 4)

#module labels as colors for nondune set
nondune_colors = labels2colors(nondune_network$colors)
#show number and size of modules
table(nondune_network$colors)

#module labels as colors for dune set
dune_colors = labels2colors(dune_network$colors)
#matching labels to those used in nondune set 
dune_colors = matchLabels(dune_colors, nondune_colors)
#show number and size of modules
table(dune_network$colors)

#average module size
mean(table(nondune_network$colors))
mean(table(dune_network$colors))

#save results 
save(nondune_network, nondune_colors, dune_network, dune_colors, file = "data/Rdata/independent_networks.Rdata")

##########################################################################################
#plot dendrograms 
#nondune
pdf(file = "plots/dendrograms/nondune_dendrogram.pdf", wi=30, he=6)
plotDendroAndColors(nondune_network$dendrograms[[1]], cbind(nondune_colors, dune_colors),
                    c("Nondune colors", "Dune colors"),
                    main = "Nondune gene dendrogram and module colors",
                    dendroLabels = FALSE)
dev.off()

#dune plotted with module labels matched to nondune labels 
pdf(file = "plots/dendrograms/dune_dendrogram.pdf", wi=30, he=6)
plotDendroAndColors(dune_network$dendrograms[[1]], cbind(dune_colors, nondune_colors),
                    c("Dune colors", "Nondune colors"),
                    main = "Dune gene dendrogram and module colors",
                    dendroLabels = FALSE)
dev.off()
