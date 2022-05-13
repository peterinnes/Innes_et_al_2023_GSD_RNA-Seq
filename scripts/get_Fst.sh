#### this script calculates per site Fst as per Weir and Cockerham 1984, using vcftools ####

VCF=~/gsd_RNA-seq/data/hc_out/dune_non-dune.nHet_filtered.hc.vcf.gz
out_dir=~/gsd_RNA-seq/analysis/Fst
sample_list=~/gsd_RNA-seq/sample_list.txt

cd $out_dir
if [ ! -e dune_samples.txt ]; then
    cat $sample_list | tail -n +2 | grep -v 'non-dune' | cut -f2 > dune_samples.txt
fi

if [ ! -e non-dune_samples.txt ]; then
        cat $sample_list | grep 'non-dune' | cut -f2 > non-dune_samples.txt
fi

echo 'calculating per-site Fst'
vcftools --gzvcf $VCF --weir-fst-pop dune_samples.txt --weir-fst-pop non-dune_samples.txt --stdout \
    > Fst_per_site.weir.txt

echo 'calculating windowed Fst'
vcftools --gzvcf $VCF --weir-fst-pop dune_samples.txt --weir-fst-pop non-dune_samples.txt --fst-window-size 1000 --fst-window-step 1000 --stdout \
    > Fst_1Kb_windows.weir.txt

