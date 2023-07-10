java -jar snpEff.jar build -gtf22 -v Ha412HOv2

cd ~/gsd_RNA-seq/analysis/snpEff
snpEff Ha412HOv2 ~/gsd_RNA-seq/data/hc_out/dune_non-dune.nHet_filtered.hc.vcf.gz > dune_non-dune.nHet_filtered.snpEff.vcf

grep 'splice' dune_non-dune.nHet_filtered.snpEff.vcf | cut -f1-2 | grep -Ff - ~/gsd_RNA-seq/analysis/Fst/Fst_per_site.weir.txt > splice_variant_Fst.txt
