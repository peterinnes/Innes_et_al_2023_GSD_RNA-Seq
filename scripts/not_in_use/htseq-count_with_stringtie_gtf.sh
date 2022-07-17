# run htseq-count to quantify expression for the DESeq2 pipeline,
#this time using the gtf annotations generated with stringtie (and corresponding STAR bams) 

gff=~/gsd_RNA-seq/data/stringtie/merged_transcripts.stringtie.gtf
bams=$(ls ~/gsd_RNA-seq/data/STAR_bams/with_stringtie_gtf/*Aligned*.bam | paste -sd' ') #using sorted bams output from STAR 

cd ~/gsd_RNA-seq/analysis/DESeq2/htseq-count_out/
htseq-count -f bam -s reverse $bams $gff > htseq-count_results.with_stringtie_gtf.2022-6-26.txt
