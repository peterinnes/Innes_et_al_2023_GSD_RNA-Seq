remotes::install_github("emjosephs/quaint")
library(quaint)

#### Expression Qpc ####
# get expression data. unsure if we should use normalized counts from DESeq2, or rlog-transformed counts from DESeq2...

load("analysis/DESeq2/DESeq2_results.Robj")
normalized_counts <- as.data.frame(DESeq2::counts(DE, normalized=TRUE)) %>%
  
df1 <- data.frame(t(normalized_counts))
#df1 <- data.frame(t(scale(normalized_counts, scale = F))) #need to have individuals in rows and genes in columns. Should we mean center the expression values by using scale(scale=F)? I think the calcQpc function does this automatically.

gene_names <- names(df1)
# read in SNP data. individuals are rows, loci are columns. instead of genotypes being coded 0-1-2 (i.e. allele counts), we have them coded as frequencies (0, 0.5, 1).
G <- read.table("data/hc_out/dune_non-dune.nHet_filtered.UTR_only.hc.freqs", header = F) #SNPs from UTRs only
G <- read.table("data/hc_out/dune_non-dune.nHet_filtered.no_inversions.hc.freqs", header = F) #SNPs from non-inversion regions only.
G <- data.matrix(G)

#G <- read.table("data/hc_out/dune_non-dune.nHet_filtered.hc.freqs", header = F) #all transcriptome SNPs
random_G <- data.matrix(sample(G, 50000)) #randomly sample 50,000 transcriptome SNPs and convert to a matrix

# get the kinship matrix
k <- make_k(myG = G) #get the kinship matrix
#k <- make_k(myG = random_G) #get the kinship matrix


# get eigen values and vectors of kinship matrix
E <- eigen(k) 
E_values <- E$values
E_vectors <- E$vectors

# Testing for selection on just first PC because this is the only PC that strongly separates dune vs. non-dune. Later PCs represent within-ecotype variation.
M <- 1

# Using the rest the of PCs to estimate Va
L <- 2:dim(k)[1]


# test for selection on each gene
allGeneOutput <- sapply(1:ncol(df1), function(x){
  myQpc = calcQpc(myZ = df1[,x], myU = E_vectors, myLambdas = E_values, myL = L, myM = M)
  return(c(gene_names[x],myQpc))
})

# look at significant results
pvals <- matrix(unlist(allGeneOutput[5,]), ncol=1, byrow=T) #get pvalues. each row corresponds to a gene, columns are to PCs. We are only interested in testing for selection with the first PC.
#pvals <- matrix(unlist(allGeneOutput[5,]), ncol=5, byrow=T) # set ncol=5 if testing for selection on first 5 PCs. 

pfdr <- data.frame(apply(pvals,2, function(x){p.adjust(x, method='fdr')})) #get adjusted pvalues
names(pfdr) <- c("PC1")
#names(pfdr) <- c("PC1", "PC2", "PC3", "PC4", "PC5")

numsig <- apply(pfdr, 2, function(x){sum(x<0.1)}) #genes with fdr <.1 across the PCs

sig_genes <- data.frame(Ha412_gene=unlist(allGeneOutput[1,]), pfdr_PC1=pfdr$PC1) %>%
  filter(pfdr_PC1 < 0.1) %>%
  mutate(Ha412_gene=gsub("mRNA.","gene:",Ha412_gene)) %>%
  left_join(Ha412_Ath_mappings)
  

write.table(sig_genes, "analysis/Qpc/sig_genes_PC1.txt", sep = '\t', row.names = F, quote = F)

# plot expression vs PC score for specific gene
# read in population data
# read in expression data

grep("mRNA.Ha412HOChr11g0499171", names(df1)) #find index of the gene we want to plot. 18680

qpc_expr_and_pops <- data.frame(pop=pop_code, gene_expr=df1[,18680])
qpc_expr_and_pops$gene_expr <- qpc_expr_and_pops$gene_expr[1:24] - mean(qpc_expr_and_pops$gene_expr) #mean center the expression values

# get eigenvectors, merge, and get eigen values
qpc_PC1 <- data.frame(PC1=E_vectors[,1], stringsAsFactors = F)
#qpc_PC2 <- data.frame(PC2=E_vectors[,2], stringsAsFactors = F)
#qpc_PC3 <- data.frame(PC3=E_vectors[,3], stringsAsFactors = F)
#qpc_PC4 <- data.frame(PC4=E_vectors[,4], stringsAsFactors = F)
#qpc_PC5 <- data.frame(PC5=E_vectors[,5], stringsAsFactors = F)


qpc_expr_and_pops <- dplyr::bind_cols(qpc_expr_and_pops[-24,], qpc_PC1) #remove last row (last individual) because this occurs in construction of kinship matrix k.
#qpc_expr_and_pops$gene_expr <- qpc_expr_and_pops$gene_expr[1:dim(E_vectors)[1]] - mean(qpc_expr_and_pops$gene_expr) #mean center the expression values


lambda <- E_values[1] #get lambda for PC1

# calculate the CI
gene_results <- allGeneOutput[,18680]
gene_name <- names(df1)[18680]
Va_est <- var0(gene_results$cml)
CI <- 1.96*sqrt(Va_est*lambda)

# linear regression of expression ~ pc1
qpc_lm <- lm(qpc_expr_and_pops$gene_expr ~ qpc_expr_and_pops$PC1)
coeff <- qpc_lm$coefficients[[2]]

plot_qpc <- ggplot(data=qpc_expr_and_pops, aes(x=PC1, y=gene_expr, color=pop, shape=pop)) +
  scale_color_manual(values=c("gold2", "forestgreen")) +
  scale_shape_manual(values=c(17,19)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        text = element_text(size=15),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=18),
        axis.title.x  = element_text(size=18),
        legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 20),
        plot.title = element_text(hjust = 0.5)) +
  geom_point(size=4, alpha=0.5) +
  geom_abline(slope = CI, linetype = 2, col = 'grey50', size=1) +
  geom_abline(slope = coeff, size = 1.25) +
  geom_abline(slope = -CI, linetype = 2, col = 'grey50', size=1)

#### Splicing Qpc ####
all_AS_events_PSI <- read.table("analysis/rMATS/results_2022-07-14/all_AS_events_PSI.txt", header = T)
rownames(all_AS_events_PSI) <- paste0(all_AS_events_PSI$GeneID,"_",all_AS_events_PSI$ID) 
df2 <- data.frame(t(all_AS_events_PSI[,-c(1:3)]))
df2 <- df2[ , apply(df2, 2, function(x) !any(is.na(x)))] #remove all events with at least one NA sample
AS_event_names <- names(df2)


allEventOutput <- sapply(1:ncol(df2), function(x){
  myQpc = calcQpc(myZ = df2[,x], myU = E_vectors, myLambdas = E_values, myL = L, myM = M)
  return(c(AS_event_names[x],myQpc))
})

AS_Qpc_pvals <- matrix(unlist(allEventOutput[5,]), ncol=1, byrow=T) #get pvalues. each row corresponds to a gene, columns are to PCs. We are only interested in testing for selection with the first PC.
#pvals <- matrix(unlist(allGeneOutput[5,]), ncol=5, byrow=T) # set ncol=5 if testing for selection on first 5 PCs. 

AS_Qpc_pfdr <- data.frame(apply(AS_Qpc_pvals,2, function(x){p.adjust(x, method='fdr')})) #get adjusted pvalues
names(AS_Qpc_pfdr) <- c("PC1")
#names(pfdr) <- c("PC1", "PC2", "PC3", "PC4", "PC5")

AS_Qpc_numsig <- apply(AS_Qpc_pfdr, 2, function(x){sum(x<0.1)}) #genes with fdr <.1 across the PCs. Zero splice events under selection...

# Using leafcutter intron excision ratios as the quantitative trait
df3 <- ier_df[-c(1:4)]

allIntronOutput <- sapply(1:ncol(df3), function(x){
  myQpc = calcQpc(myZ = df3[,x], myU = E_vectors, myLambdas = E_values, myL = L, myM = M)
  return(c(AS_event_names[x],myQpc))
})

leafcutter_Qpc_pvals <- matrix(unlist(allIntronOutput[5,]), ncol=1, byrow=T) #get pvalues. each row corresponds to a gene, columns are to PCs. We are only interested in testing for selection with the first PC.
#pvals <- matrix(unlist(allGeneOutput[5,]), ncol=5, byrow=T) # set ncol=5 if testing for selection on first 5 PCs. 

leafcutter_Qpc_pfdr <- data.frame(apply(leafcutter_Qpc_pvals,2, function(x){p.adjust(x, method='fdr')})) #get adjusted pvalues
names(leafcutter_Qpc_pfdr) <- c("PC1")
#names(pfdr) <- c("PC1", "PC2", "PC3", "PC4", "PC5")

leafcutter_Qpc_numsig <- apply(leafcutter_Qpc_pfdr, 2, function(x){sum(x<0.1)}) #genes with fdr <.1 across the PCs. Zero intron excision ratios under selection.


