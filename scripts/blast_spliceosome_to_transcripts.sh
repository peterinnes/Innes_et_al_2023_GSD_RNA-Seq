splice_components_prot=~/gsd_RNA-seq/analysis/spliceosome_prot_components_combined.aa.fasta
splice_components_ncRNA=~/gsd_RNA-seq/analysis/spliceosome_ncRNA_components_combined.nucl.fasta
prot_blastdb=~/gsd_RNA-seq/data/ref_genome_Ha412HO/HAN412_Eugene_curated_v1_1_cds.fa
RNA_blastdb=~/gsd_RNA-seq/data/ref_genome_Ha412HO/HAN412_Eugene_curated_v1_1_exon_merge.fa
outdir=~/gsd_RNA-seq/analysis/BLAST_out

tblastn -query $splice_components_prot -db $prot_blastdb -num_threads 12 -outfmt 6 -out $outdir/spliceosome_components_prot_vs_HAN412.tblastn
blastn -query $splice_components_ncRNA -db $RNA_blastdb -num_threads 12 -outfmt 6 -out $outdir/spliceosome_components_ncRNA_vs_HAN412.blastn
