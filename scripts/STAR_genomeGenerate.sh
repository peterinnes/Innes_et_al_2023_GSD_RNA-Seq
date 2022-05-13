# Generate genome index file for STAR

STAR \
	--runThreadN 16 \
	--runMode genomeGenerate \
	--genomeDir /data3/peter/sunflower/ref_genome_Ha412HO/STAR_genome_indices \
	--genomeFastaFiles /data3/peter/sunflower/ref_genome_Ha412HO/Ha412HOv2.0-20181130.fasta \
	--sjdbGTFfile /data3/peter/sunflower/ref_genome_Ha412HO/HAN412_Eugene_curated_v1_1.gff3 \
	--sjdbGTFtagExonParentTranscript Parent \
	--sjdbOverhang 74 
