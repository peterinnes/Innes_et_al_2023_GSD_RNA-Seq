rmats.py --b1 ~/gsd_RNA-seq/data/STAR_bams/dune_bams_list.txt --b2 ~/gsd_RNA-seq/data/STAR_bams/non-dune_bams_list.txt
 --gtf ~/gsd_RNA-seq/data/ref_genome_Ha412HO/HAN412_Eugene_curated_v1_1.gtf -t single --libType fr-firststrand --readLength 75 --variable-read-length --allow-clipping --n
ovelSS --cstat .05 --nthread 1 --od results_2022-05-24/ --tmp results_2022-05-24/tmp/
