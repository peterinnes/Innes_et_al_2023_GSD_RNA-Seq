fqdir=~/gsd_RNA-seq/data/fastq/trimmed/

cd $fqdir
for fq in *fq*; do
	sample=$(echo $fq | cut -d '_' -f1) 
	echo "quantifying $sample..."
	kallisto quant \
		-i ~/gsd_RNA-seq/data/transcriptome/multi_sample_Trinity.kallisto.idx \
		-t 8 \
		--rf-stranded \
		-o ~/gsd_RNA-seq/data/kallisto_out/$sample \
		--single \
		-l 180 \
		-s 20 \
		$fq
done
