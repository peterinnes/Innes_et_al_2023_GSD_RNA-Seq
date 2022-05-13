# Mapping RNA-seq reads to Sunflower reference GENOME using STAR

cd ~/gsd_RNA-seq/data/fastq/trimmed/
for fq in $(ls Kane-602-*); do
	sample=$(echo $fq | cut -d '_' -f1-2)
	STAR \
		--runThreadN 16 \
		--twopassMode Basic \
		--genomeDir /data3/peter/sunflower/ref_genome_Ha412HO/STAR_genome_indices \
		--readFilesIn $fq \
		--readFilesCommand zcat \
		--outFileNamePrefix ~/gsd_RNA-seq/data/STAR_bams/${sample}_ \
		--outSAMtype BAM SortedByCoordinate
done
