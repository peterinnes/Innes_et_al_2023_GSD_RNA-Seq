#WGCNA (weighted gene co-expression network analysis) for Kane lab rotation project
#By Kaylee Rosenberger, following WGCNA package tutorial 
#using Great Sand Dunes National Park prairie sunflower genome (from Peter Innes)
#dataset is a table of 30000+ genes for 25 samples of dune and non-dune populations 
#each cell represents regularized log transformed expression counts 

#This script compares intramodular connectivity for DEG and nonDEG

#import libraries
library(BiocManager)
library(WGCNA)
library(tibble)
library(dplyr)
library(ggplot2)

#############################################################################################
#load in connectivity info
load("data/Rdata/whole_network_connectivity.Rdata")

#############################################################################################
#calculating change in connectivity from nondune to dune 

dune_WN_con = dune_WN_con[,c(2,5,7,8)]
colnames(dune_WN_con) = c("IM_dune", "mod_dune", "gene", "diff_exp")
non_dune_WN_con = non_dune_WN_con[,c(2,5)]
colnames(non_dune_WN_con) = c("IM_nondune", "mod_nondune")

conn_change = cbind(non_dune_WN_con, dune_WN_con)
#removing genes not assigned to moduels in dune and nondune
conn_change = subset(conn_change, !is.na(IM_nondune))
conn_change = subset(conn_change, !is.na(IM_dune))

#negative means more connected in dune
#positive means more connected in nondune 
conn_change$change = conn_change$IM_nondune - conn_change$IM_dune

#plotting results with lines indicating the top 10% genes with changing connectivity
#in both directions (positive and negative)
ggplot(conn_change, aes(x=change)) +
  geom_histogram(binwidth = 5, color = "black", fill = "lightblue") +
  xlab("Change in intramodular connectivity") +
  ylab("Frequency") +
  ggtitle("Change in connectivty from nondune to dune network") +
  theme_bw() +
  geom_vline(aes(xintercept = quantile(conn_change$change,0.95)), linetype = "dashed") +
  geom_vline(aes(xintercept = quantile(conn_change$change,0.05)), linetype = "dashed")

#getting top 10% of genes changing connectivity 
conn_change_pos = conn_change %>% filter(quantile(change, 0.95)<change)
conn_change_neg = conn_change %>% filter(quantile(change,0.05)>change)

conn_change_top10 = rbind(conn_change_pos, conn_change_neg)

save(conn_change, conn_change_neg, conn_change_pos, conn_change_top10, file = "data/Rdata/change_connectivity.Rdata")
