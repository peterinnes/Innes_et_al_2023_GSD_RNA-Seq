# script to make gvcf files for each sample and then joint-call genotypes for a single, final vcf table (unfiltered)
# based on GATK RNA-seq best practices: https://github.com/gatk-workflows/gatk4-rnaseq-germline-snps-indels

ref=~/gsd_RNA-seq/data/ref_genome_Ha412HO/Ha412HOv2.0-20181130.fasta
bam_dir=~/gsd_RNA-seq/data/STAR_bams/
out_dir=~/gsd_RNA-seq/data/hc_out/

ls ${bam_dir}*sorted_RGadded_dupmarked_split.bam | while read bam; do
    sample=`echo $bam | cut -d "/" -f7 | cut -d "." -f1`
    echo "Current file is $bam, sample is $sample..."

    # make gvcf file. chase confidence threshold per GATK RNA-seq best practices
    gatk HaplotypeCaller \
        -R $ref \
        -I $bam \
        -O $out_dir$sample.g.vcf.gz \
        --emit-ref-confidence GVCF \
	--dont-use-soft-clipped-bases \
	--standard-min-confidence-threshold-for-calling 20
done

# get list of gvcfs for input string
for gvcf in $( ls $out_dir*.g.vcf.gz ); do
    gvcf_args+="--variant $gvcf "
done

echo "finished making gvcfs: $gvcf_args"

# combine gvcfs
echo "combining gvcfs...$gvcf_args"
gatk CombineGVCFs \
    -R $ref \
    $gvcf_args \
    -O ${out_dir}dune_non-dune.joint.vcf.gz

# genotype!
echo "beginning final genotyping..."
gatk GenotypeGVCFs \
    -R $ref \
    -V ${out_dir}dune_non-dune.joint.vcf.gz \
    -O ${out_dir}dune_non-dune.vcf.gz 

# clean up
#rm $out_dir*.g.vcf*
#rm $out_dir*.joint.vcf.gz
