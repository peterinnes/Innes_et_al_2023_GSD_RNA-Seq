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

echo 'calculating 1kb windowed Fst'
vcftools --gzvcf $VCF --weir-fst-pop dune_samples.txt --weir-fst-pop non-dune_samples.txt --fst-window-size 1000 --fst-window-step 1000 --stdout \
    > Fst_1Kb_windows.weir.txt

echo 'calculating 10kb windowed Fst'
vcftools --gzvcf $VCF --weir-fst-pop dune_samples.txt --weir-fst-pop non-dune_samples.txt --fst-window-size 10000 --fst-window-step 10000 --stdout \
    > Fst_10Kb_windows.weir.txt

echo 'calculating 100kb windowed Fst'
vcftools --gzvcf $VCF --weir-fst-pop dune_samples.txt --weir-fst-pop non-dune_samples.txt --fst-window-size 100000 --fst-window-step 100000 --stdout \
	> Fst_100Kb_windows.weir.txt 

echo 'calculating 1Mb windowed Fst'
vcftools --gzvcf $VCF --weir-fst-pop dune_samples.txt --weir-fst-pop non-dune_samples.txt --fst-window-size 1000000 --fst-window-step 1000000 --stdout \
    > Fst_1Mb_windows.weir.txt
echo 'calculating 2Mb windowed Fst'
vcftools --gzvcf $VCF --weir-fst-pop dune_samples.txt --weir-fst-pop non-dune_samples.txt --fst-window-size 2000000 --fst-window-step 2000000 --stdout \
    > Fst_2Mb_windows.weir.txt
