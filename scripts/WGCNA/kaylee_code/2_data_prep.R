#WGCNA (weighted gene co-expression network analysis) for Kane lab rotation project
#By Kaylee Rosenberger, following WGCNA package tutorial 
#using Great Sand Dunes National Park prairie sunflower genome (from Peter Innes)
#dataset is a table of 30000+ genes for 25 samples of dune and non-dune populations 
#each cell represents regularized log transformed expression counts 

#This script reformats expression data to separate dune and non-dune ecotypes 
#for a consensus analysis
#It also checks for genes with too many missing samples, which should be removed
#lastly it checks the clustering of samples to determine if there are any outliers
#that could impact network construction 

#import libraries
library(BiocManager)
library(WGCNA)
library(dplyr)

############################################################################################
#Load in data 
options(stringsAsFactors = FALSE) #necessary setting 
expression_data = read.csv("data/expData/rlog-transformed_expression_counts.tsv", sep="\t")
# View(expression_data) #view data to check import 
# dim(expression_data)
# names(expression_data)

############################################################################################
#Filter data

#filter data for dune vs. non-dune types to create two separate dataframes to work with 
non_dune_expression_data = expression_data %>% select(starts_with(c("Ha412", "non.dune")))
dune_expression_data = expression_data %>% select(starts_with(c("Ha412", "dune")))

#converting data to multi-set format suitable for consensus analysis
num_sets = 2 #two sets--dune and non dune
set_labels = c("Non-dune", "Dune") #for labeling plots 
multi_exp = vector(mode = "list", length = num_sets) #creating a vector of each dataset
#working with non dune data first
#transpose data and remove first column, which has the genes 
multi_exp[[1]] = list(data = as.data.frame(t(non_dune_expression_data[,2:13]))) 
colnames(multi_exp[[1]]$data) = non_dune_expression_data$Ha412_gene #getting gene names for the columns
rownames(multi_exp[[1]]$data) = names(non_dune_expression_data[,2:13]) #getting sample names for the rows
#working with non dune data first
#transpose data and remove first column, which has the genes listed
multi_exp[[2]] = list(data = as.data.frame(t(dune_expression_data[,2:13])))
colnames(multi_exp[[2]]$data) = dune_expression_data$Ha412_gene
rownames(multi_exp[[2]]$data) = names(dune_expression_data[,2:13])
#check to make sure the data is formatted properly 
expr_size = checkSets(multi_exp)

###########################################################################################
#Check for genes with too many missing samples
gsg = goodSamplesGenesMS(multi_exp, verbose = 3)
gsg$allOK #if TRUE all genes pass the cut 

###########################################################################################
#Check clustering of samples to determine if there are any outliers 
#cluster samples based on Euclidean distance 
sampleTrees = list()
for (set in 1:num_sets) {
  sampleTrees[[set]] = hclust(dist(multi_exp[[set]]$data), method = "average")
}
#plotting cluster dendrograms of samples 
pdf(file = "../plots/sampleClustering.pdf", width = 12, height = 12);
par(mfrow=c(2,1))
par(mar = c(0, 4, 2, 0))
for(set in 1:num_sets) {
  plot(sampleTrees[[set]], main = paste("Sample clustering on all genes in", set_labels[set]),
       xlab="", sub="", cex = 0.7)
}
dev.off()
#check plots for obvious outliers, and cut the tree if necessary 
#But it's not necessary here

###########################################################################################
#save data
save(multi_exp, expr_size, set_labels, file = "data/Rdata/multi_expr.Rdata")
