# find specific matches between Trinity 'genes' (just using longest isoform of each 'gene') and annotated HA412HOV2 genes.

blastn -query ~/gsd_RNA-seq/data2/transcriptome/all_sample_Trinity.cd-hit-est_99.longest_isos.single_line.fasta -db ~/gsd_RNA-seq/data/ref_genome_Ha412HO/HAN412_Eugene_curated_v1_1_all_transcripts.fa -outfmt 6 > ~/gsd_RNA-seq/analysis/BLAST_out/Trinity_transcripts_vs_HAN412_transcripts.v2.blast

