ref=~/gsd_RNA-seq/data/ref_genome_Ha412HO/Ha412HOv2.0-20181130.fasta
bam_dir=~/gsd_RNA-seq/data/STAR_bams/
outdir=~/gsd_RNA-seq/data/mpileup_out/

ls $bam_dir*_Aligned*.bam > ${bam_dir}bam_list.txt

bcftools mpileup \
    -Ou \
    -a DP \
    -f $ref \
    -b ${bam_dir}bam_list.txt | \
    #-q 20 | \
bcftools call \
    -Oz \
    -f GQ \
    --variants-only \
    --multiallelic-caller > ${outdir}dune_non-dune_2022-2-5.unfiltered.mpileup.vcf.gz
