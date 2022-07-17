## Differential Alternative Splicing analysis using ASPli. NOTE: ended up not pursuing this analysis because #ASPlie requires splicing to be annotated in the gff file (i.e. multiple isoforms annotated per gene). #Unfortunately the Ha412v2 assembly only has one isoform annotated per gene.
##
## vignette(ASPli)
#
library(ASpli)
library(GenomicFeatures)
#
## get input data
TxDb <- makeTxDbFromGFF(file="data/stringtie/merged_transcripts.stringtie.gtf", format = c("gtf"),organism = "Helianthus petiolaris")
features <- binGenome(TxDb)
geneCoord <- featuresg( features )
binCoord <- featuresb( features )
junctionCoord <- featuresj( features )
                              
bam_files <- read.table("data/STAR_bams/with_stringtie_gtf/bam_list.txt", col.names = "bam")
samples <- read.table("samples.txt", col.names = "sample")
targets <- data.frame(bam=bam_files, habitat=sample_table_deseq2$habitat)
rownames(targets)<- samples$sample
getConditions(targets)
#
# gbCounts: Summarize read overlaps against all feature levels i.e. Read counting against annotated features
gbcounts <- gbCounts( features = features,
                      targets = targets,
                      minReadLength = 75, maxISize = 50000,
                      libType="SE",
                      strandMode = 2)

# jCounts: Summarize junctions inclusion indices PSI, PIR and PJU i.e. 
# Junction-based de-novo counting and splicing signal estimation:
jcounts <- jCounts(counts = gbcounts,
                 features = features,
                 minReadLength = 75,
                 libType="SE",
                 strandMode=2)

# Bin-based coverage differential signals of AS. Can filter for low-expressed genes with minGenReads and minBinReads. I should probably match the same filtering as used w/ rMATS: 24 / 12 reads.
gb <- gbDUreport(gbcounts, contrast = c(1,-1),
                  minGenReads = 24,
                  minBinReads = 12)
                  #minRds = 0.05,
                  #ignoreExternal = TRUE,
                  #ignoreIo = TRUE,
                  #ignoreI = FALSE,
                  #filterWithContrasted = TRUE,
                  #verbose = TRUE,
                  #formula = NULL,
                  #coef = NULL)
geneX <- genesDE( gb )
binDU <- binsDU( gb )

writeDU(gb, output.dir = "analysis/ASpli") #export results as tab-delimited table

# Junction-centered analysis of AS
jdu <- jDUreport(jcounts, contrast = c(1,-1),
                 minAvgCounts = 12)
                 #filterWithContrasted = TRUE,
                 #runUniformityTest = FALSE,
                 #maxPValForUniformityCheck = 0.2,
                 #strongFilter = TRUE,
                 #maxConditionsForDispersionEstimate = 24,
                 #formula = NULL,
                 #coef = NULL,
                 #maxFDRForParticipation = 0.05,
                 #useSubset = FALSE)
dim(localej( jdu ))
dim(localec( jdu ))
dim(anchorj( jdu ))
dim(anchorc( jdu ))
dim(jir( jdu ))
dim(jes( jdu ))
dim(jalt( jdu ))

# bin and junction signal integration
sr <- splicingReport(gb, jdu, counts=gbcounts)
writeSplicingReport(sr, output.dir = "analysis/ASpli")

integrate_splice_signals <- integrateSignals(sr, jcounts)
browseVignettes("ASpli")
