# filtering our variants

VCF_IN=~/gsd_RNA-seq/data/mpileup_out/dune_non-dune_2022-2-5.unfiltered.mpileup.vcf.gz
VCF_OUT=~/gsd_RNA-seq/data/mpileup_out/dune_non-dune_2022-2-5.filtered.mpileup.vcf.gz

# set filters
MAF=0.05
MISS=.75
QUAL=30
MIN_DEPTH=3
MAX_DEPTH=20 #probably shouldn't use max depth filter for RNA-seq data

vcftools --gzvcf $VCF_IN \
--remove-indels --min-alleles 2 --max-alleles 2 --maf $MAF --max-missing $MISS --minQ $QUAL \
--min-meanDP $MIN_DEPTH --max-meanDP $MAX_DEPTH \
--minDP $MIN_DEPTH --maxDP $MAX_DEPTH --recode --stdout | gzip -c > \
$VCF_OUT

#minDP/maxDP: min or maxdepth allowed for a genotype - any individual failing this threshold is marked as having a missing genotype. min-meanDP and max-meanDP are the min/max depth for a site.

