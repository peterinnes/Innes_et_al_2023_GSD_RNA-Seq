# count number of heterozygous indivuduals, per site. From Kevin Blighe https://www.biostars.org/p/291147/

#paste <(bcftools view data/hc_out/dune_non-dune.filtered.hc.vcf.gz |\
#    awk -F"\t" 'BEGIN {print "CHR\tPOS\tID\tREF\tALT"} \
#      !/^#/ {print $1"\t"$2"\t"$3"\t"$4"\t"$5}') \
#    \
#  <(bcftools query -f '[\t%SAMPLE=%GT]\n' data/hc_out/dune_non-dune.filtered.hc.vcf.gz |\
#    awk 'BEGIN {print "nHet"} {print gsub(/0\|1|1\|0|0\/1|1\/0/, "")}')

# alternatively, use vcftools --hardy (it outputs heterozygosity stats per site in addition to testing for HWE)
vcftools --gzvcf data/hc_out/dune_non-dune.filtered.hc.vcf.gz --hardy --output hwe_per_site.txt
