# Mapping RNA-seq reads to Sunflower reference GENOME using STAR, this using the gtf file we generated with string tie, which (hopefully) has more complete annotation of additional isoforms per gene.

cd ~/gsd_RNA-seq/data/fastq/trimmed/
for fq in $(ls Kane-602-*); do
	sample=$(echo $fq | cut -d '_' -f1-2)
	STAR \
		--runThreadN 8 \
		--twopassMode Basic \
		--genomeDir /data3/peter/sunflower/stringtie/STAR_genome_indices \
		--readFilesIn $fq \
		--readFilesCommand zcat \
		--outFileNamePrefix ~/gsd_RNA-seq/data/STAR_bams/with_stringtie_gtf/${sample}_ \
		--outSAMtype BAM SortedByCoordinate
done
