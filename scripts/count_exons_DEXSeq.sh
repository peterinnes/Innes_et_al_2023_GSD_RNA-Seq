gff=/home/peter/gsd_RNA-seq/data/ref_genome_Ha412HO/HAN412_Eugene_curated_v1_1.DEXSeq.gff
bams=`ls ~/gsd_RNA-seq/data/STAR_bams/*Aligned*.bam | paste -sd' '`
out_dir=~/gsd_RNA-seq/data/dexseq-count_out/

# prepare the annotation. $gff is our newly formated annotations file. 
#python /home/peter/R/x86_64-pc-linux-gnu-library/4.1/DEXSeq/python_scripts/dexseq_prepare_annotation.py HAN412_Eugene_curated_v1_1.gtf $gff

for bam in $bams; do
	sample=$(echo $bam | cut -d '/' -f7 | cut -d '_' -f1-2)
	# count
	echo "counting for $sample..."
	python3.7 /home/peter/R/x86_64-pc-linux-gnu-library/4.1/DEXSeq/python_scripts/dexseq_count.py \
		--stranded=reverse \
		--format=bam \
		$gff \
		$bam \
		$out_dir$sample.dexseq_counts.txt
done
	
