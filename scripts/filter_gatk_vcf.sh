ref=~/gsd_RNA-seq/data/ref_genome_Ha412HO/Ha412HOv2.0-20181130.fasta
vcf=~/gsd_RNA-seq/data/hc_out/dune_non-dune.unfiltered.hc.vcf.gz
out_dir=~/gsd_RNA-seq/data/hc_out/

# remove indels, select only bi-allelic SNPs
gatk SelectVariants \
 -R $ref \
 -V $vcf \
 --select-type-to-include SNP \
 --restrict-alleles-to BIALLELIC \
 -O ${out_dir}dune_non-dune.biallelic_SNPs.hc.vcf.gz

# generic hard filters for RNAseq data, per GATK:https://github.com/gatk-workflows/gatk4-rnaseq-germline-snps-indels/blob/master/gatk4-rna-best-practices.wdl
gatk VariantFiltration \
 -R $ref \
 -V ${out_dir}dune_non-dune.biallelic_SNPs.hc.vcf.gz \
 --window 35 \
 --cluster 3 \
 --filter-name "FS" \
 --filter "FS > 30.0" \
 --filter-name "QD" \
 --filter "QD < 2.0" \
 -O ${out_dir}dune_non-dune.gatk_filtered.hc.vcf.gz

# filters for vcftools
MAF=0.05
MISS=1
QUAL=30
MIN_DEPTH=5

vcftools --gzvcf ${out_dir}dune_non-dune.gatk_filtered.hc.vcf.gz \
--maf $MAF --max-missing $MISS --minQ $QUAL \
--min-meanDP $MIN_DEPTH \
--minDP $MIN_DEPTH --recode --stdout | gzip -c \
> ${out_dir}dune_non-dune.filtered.hc.vcf.gz
