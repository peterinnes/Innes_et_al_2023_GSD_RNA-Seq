# make bed file with gene positions
agat_convert_sp_gff2bed.pl \
	--gff HAN412_Eugene_curated_v1_1.gff3 \
	-sub gene \
	-o HAN412_Eugene_curated_v1_1_genes.bed

# single gene windows with no buffer
#cut -f1-3 HAN412_Eugene_curated_v1_1_genes.bed > ~/gsd_RNA-seq/analysis/Fst/single_gene_windows.bed

# single gene windows with 5kb buffer. This is the input for calculate_Fst_per_gene.py
awk 'BEGIN{OFS="\t"} {$2=$2-5000;$3=$3+5000;print $1,$2,$3}' HAN412_Eugene_curated_v1_1_genes.bed \
	> ~/gsd_RNA-seq/analysis/Fst/single_gene_windows_5kb_buffer.bed
