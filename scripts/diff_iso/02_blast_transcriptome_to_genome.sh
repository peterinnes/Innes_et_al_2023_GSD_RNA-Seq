transcriptome=~/gsd_RNA-seq/data2/transcriptome/all_sample_Trinity.cd-hit-est_99.fasta
genome=~/gsd_RNA-seq/data/ref_genome_Ha412HO/Ha412HOv2.0-20181130.fasta
outdir=~/gsd_RNA-seq/analysis/BLAST_out

#makeblastdb -in $genome -dbtype nucl
blastn -query $transcriptome -db $genome -num_threads 12 -outfmt "6 qseqid sseqid pident qlen length qstart qend sstart send evalue gaps" -perc_identity 85 -out $outdir/petiolaris_transcriptome_vs_annuus_genome.blast #don't go lower than 85
