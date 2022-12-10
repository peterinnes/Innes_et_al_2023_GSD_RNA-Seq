##########################################################################################
#import libraries
library(BiocManager)
library(WGCNA)
library(tibble)
library(dplyr)

#load in data
load("data/Rdata/independent_networks.Rdata")

##########################################################################################
#running module preservation analysis 
#comparing nondune to dune 
#how well are the nondune modules preserved in the dune datasets? 
multi_color = list("Non-dune" = nondune_colors, "Dune" = dune_colors)

module_preservation = modulePreservation(multi_exp, multi_color, dataIsExpr = TRUE, 
                                         networkType = "signed", referenceNetworks = c(1:2),
                                         nPermutations = 200, randomSeed = 1,
                                         quickCor = 0, verbose = 3)
#save file because it takes a long time to run 
save(module_preservation, file = "module_preservation_stats.Rdata")

#########################################################################################
#Analysis and visualization of module preservation 
load("data/Rdata/module_preservation_stats.Rdata")

ref = 2
test = 1
statsObs = cbind(module_preservation$quality$observed[[ref]][[test]][, -1], module_preservation$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(module_preservation$quality$Z[[ref]][[test]][, -1], module_preservation$preservation$Z[[ref]][[test]][, -1])

print( cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
             signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) )

#plotting preservation
modColors = rownames(module_preservation$preservation$observed[[ref]][[test]])
moduleSizes = module_preservation$preservation$Z[[ref]][[test]][, 1];
# leave grey and gold modules out
plotMods = !(modColors %in% c("grey", "gold"));
# Text labels for points
text = modColors[plotMods];
# Auxiliary convenience variable
plotData = cbind(module_preservation$preservation$observed[[ref]][[test]][, 2], module_preservation$preservation$Z[[ref]][[test]][, 2])
# Main titles for the plot
mains = c("Preservation Median rank", "Preservation Zsummary");
# Start the plot
pdf(file="plots/module_preservation.pdf", wi=10, h=5)
par(mfrow = c(1,2))
par(mar = c(4.5,4.5,2.5,1))
for(p in 1:2) {
  min = min(plotData[, p], na.rm = TRUE);
  max = max(plotData[, p], na.rm = TRUE);
  # Adjust ploting ranges appropriately
  if (p==2)
  {
    if (min > -max/10) min = -max/10
    ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
  } else
    ylim = c(max + 0.1 * (max-min), min - 0.1 * (max-min))
  plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21,
       main = mains[p],
       cex = 2.4,
       ylab = mains[p], xlab = "Module size", log = "x",
       ylim = ylim,
       xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4)
  #labelPoints(moduleSizes[plotMods], plotData[plotMods, p], text, cex = 1, offs = 0.08)
  # For Zsummary, add threshold lines
  if (p==2)
  {
    abline(h=0)
    abline(h=2, col = "blue", lty = 2)
    abline(h=10, col = "darkgreen", lty = 2)
  }
}
# If plotting into a file, close it
dev.off()

##########################################################################################
#Determining the number of moduels with the same labels 

length(intersect(unique(dune_colors), unique(nondune_colors)))

#determining unique modules
outersect <- function(x, y) {
  sort(c(setdiff(x, y),
         setdiff(y, x)))
}
outersect(unique(dune_colors), unique(nondune_colors))

#Determining unique modules to each network: 
#determining modules unique to nondune
nondune_unique_mods = setdiff(unique(nondune_colors), unique(dune_colors))

#determining modules unique to dune 
dune_unique_mods = setdiff(unique(dune_colors), unique(nondune_colors))

