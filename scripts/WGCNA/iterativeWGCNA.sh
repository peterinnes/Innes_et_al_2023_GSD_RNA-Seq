# conda activate EDTA

iterativeWGCNA -i ~/gsd_RNA-seq/analysis/WGCNA/iterativeWGCNA/dune_expr_data.iterWGCNA.tsv \
	-o ~/gsd_RNA-seq/analysis/WGCNA/iterativeWGCNA/dune \
	--verbose \
	--wgcnaParameters maxBlockSize=25000,power=18,minKMEtoStay=0.8,minCoreKME=0.8,reassignThreshold=0.05,minModuleSize=30,nThreads=20 \
	--enableWGCNAThreads

iterativeWGCNA -i ~/gsd_RNA-seq/analysis/WGCNA/iterativeWGCNA/non_dune_expr_data.iterWGCNA.tsv \
	-o ~/gsd_RNA-seq/analysis/WGCNA/iterativeWGCNA/non_dune \
	--verbose \
	--wgcnaParameters maxBlockSize=25000,power=18,minKMEtoStay=0.8,minCoreKME=0.8,reassignThreshold=0.05,minModuleSize=30,nThreads=20 \
	--enableWGCNAThreads
