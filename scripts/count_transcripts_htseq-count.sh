# run htseq-count to quantify expression for the DESeq2 pipeline. Need to set --idattr to Parent so that reads are counted for genes rather than individual exons.

bams=$(ls ~/gsd_RNA-seq/data/STAR_bams/*Aligned*.bam | paste -sd' ') #using sorted bams output from STAR
gff=~/gsd_RNA-seq/data/ref_genome_Ha412HO/HAN412_Eugene_curated_v1_1.gff3

cd ~/gsd_RNA-seq/analysis/DESeq2/htseq-count_out/
htseq-count -f bam -s reverse --idattr Parent $bams $gff > htseq-count_results.2022-6-26.txt
