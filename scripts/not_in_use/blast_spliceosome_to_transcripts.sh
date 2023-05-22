query=~/gsd_RNA-seq/spliceosome_components_combined.aa.fasta
blastdb=~/gsd_RNA-seq/data/ref_genome_Ha412HO/HAN412_Eugene_curated_v1_1_cds.fa

tblastn -query $query -db $blastdb -num_threads 12 -outfmt 6 -out ~/gsd_RNA-seq/analysis/BLAST_out/spliceosome_components_vs_HAN412.tblastn
