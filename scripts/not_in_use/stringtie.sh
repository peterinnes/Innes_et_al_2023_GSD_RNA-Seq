rundir=~/gsd_RNA-seq/data/STAR_bams/for_stringtie/
outdir=~/gsd_RNA-seq/data/stringtie/
ref=~/gsd_RNA-seq/data/ref_genome_Ha412HO/HAN412_Eugene_curated_v1_1.gff3

# assemble transcripts for each sample, using reference annotations
cd $rundir
for bam in *bam; do
	sample=$(echo $bam | cut -d "_" -f1-2)
	stringtie $bam \
		-p 8 \
		-v \
		-o $outdir$sample.stringtie.gtf \
		-G $ref \
		--rf \
		--conservative
done

# merge all the gtf files
cd $outdir
stringtie --merge \
	-G $ref \
	-o merged_transcripts.stringtie.gtf \
	gtf_list.txt
	
