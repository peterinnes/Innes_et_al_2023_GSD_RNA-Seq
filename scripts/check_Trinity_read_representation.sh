# This script aligns reads to the assembled transcriptome the same reads were used to generate, as a first step in quality assessment of the assembly
# as described: https://github.com/trinityrnaseq/trinityrnaseq/wiki/RNA-Seq-Read-Representation-by-Trinity-Assembly

assembly=
bowtie2-build $assembly $assembly

bowtie2 -p 10 -q --no-unal -k 20 -x $assembly -U $fq \ 
	2>align_stats.txt | samtools view -@10 -Sb -o $sample.bam

# then visualize the statistics:
#cat 2>&1 align_stats.txt
