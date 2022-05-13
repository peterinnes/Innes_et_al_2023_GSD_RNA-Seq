fqdir=~/data/sunflower/fastq/

for fq in $(ls ${fqdir}raw/); do
        sample=$(echo $fq | cut -d '.' -f1); echo trimming $sample...
	fastp -i ${fqdir}raw/$fq -o ${fqdir}trimmed/$sample.trimmed.fq.gz \
		-h fastp_$sample.html -j fastp_$sample.json \
		--cut_tail --cut_tail_window_size 1 \
		--dont_eval_duplication \
		--overrepresentation_analysis
done
