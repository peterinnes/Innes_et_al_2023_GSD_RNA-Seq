query=~/gsd_RNA-seq/spliceosome_components.aa.fasta
blastdb=~/gsd_RNA-seq/data/ref_genome_Ha412HO/HAN412_Eugene_curated_v1_1_cds.fa

tblastn -query $query -db $blastdb -evalue .0001 -num_threads 4 -max_target_seqs 4 -outfmt "6 qseqid sseqid pident qlen length qstart qend sstart send evalue bitscore" -out ~/gsd_RNA-seq/analysis/BLAST_out/spliceosome_components_vs_HAN412.tblastn
