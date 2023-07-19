# count number of heterozygous indivuduals, per site,
# with  vcftools --hardy 
# (it outputs heterozygosity stats per site in addition to testing for HWE)
vcftools --gzvcf data/hc_out/dune_non-dune.filtered.hc.vcf.gz \
	--hardy \
	--output nHet_out.txt

# then find sites with > 60% het. THis filters 6660 variants. 50% (>12) would filter 18907; >13 (54%) would filter 11424
awk '{if($6 > 14) {print $1"\t"$2"\t"$3"\t"$4"\t"$5}}' nHet_out.txt \
	> excess_het_loci.txt
zgrep -vFf excess_het_loci.txt dune_non-dune.filtered.hc.vcf.gz \
	> dune_non-dune.nHet_filtered.hc.vcf

# old code to manually count hets per site, 
# from Kevin Blighe https://www.biostars.org/p/291147/
#paste <(bcftools view data/hc_out/dune_non-dune.filtered.hc.vcf.gz |\
#    awk -F"\t" 'BEGIN {print "CHR\tPOS\tID\tREF\tALT"} \
#      !/^#/ {print $1"\t"$2"\t"$3"\t"$4"\t"$5}') \
#    \
#  <(bcftools query -f '[\t%SAMPLE=%GT]\n' data/hc_out/dune_non-dune.filtered.hc.vcf.gz |\
#    awk 'BEGIN {print "nHet"} {print gsub(/0\|1|1\|0|0\/1|1\/0/, "")}')


