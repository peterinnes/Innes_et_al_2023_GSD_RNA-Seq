# Mapping RNA-seq reads to Sunflower reference GENOME using STAR, this time for compatibility with stringtie, which we subsequently use to make new reference-based annotations of transcripts. Add the flag "--outSAMstrandField intronMotif" for stringtie compatibility. Per STAR manual, this option "changes the output alignments: reads with inconsistent and/or non-canonical introns are filtered out.

cd ~/gsd_RNA-seq/data/fastq/trimmed/
for fq in $(ls Kane-602-*); do
	sample=$(echo $fq | cut -d '_' -f1-2)
	STAR \
		--runThreadN 8 \
		--twopassMode Basic \
		--genomeDir /data3/peter/sunflower/ref_genome_Ha412HO/STAR_genome_indices \
		--readFilesIn $fq \
		--readFilesCommand zcat \
		--outSAMstrandField intronMotif \
		--outFileNamePrefix ~/gsd_RNA-seq/data/STAR_bams/for_stringtie/${sample}_ \
		--outSAMtype BAM SortedByCoordinate
done
