# Generate genome index file for STAR alignments, this time using the annotations produced by stringtie

STAR \
	--runThreadN 16 \
	--runMode genomeGenerate \
	--genomeDir /data3/peter/sunflower/stringtie/STAR_genome_indices \
	--genomeFastaFiles /data3/peter/sunflower/ref_genome_Ha412HO/Ha412HOv2.0-20181130.fasta \
	--sjdbGTFfile /data3/peter/sunflower/stringtie/merged_transcripts.stringtie.gtf \
	--sjdbOverhang 74 
