#### this script finds matches between Trinity 'genes' (just using longest isoform of each 'gene') and annotated HA412HOV2 genes ####
cd ~/gsd_RNA-seq/data/ref_genome_Ha412HO/

# conda activate agat

# extract Ha412HO transcript sequences with introns. Use this bc Trinity isoforms could have intronic sequence in the case of intron retention
#agat_sp_extract_sequences.pl -g HAN412_Eugene_curated_v1_1.gff3 -f Ha412HOv2.0-20181130.fasta -t rRNA -o HAN412_Eugene_curated_v1_1_rRNA.fa
#agat_sp_extract_sequences.pl -g HAN412_Eugene_curated_v1_1.gff3 -f Ha412HOv2.0-20181130.fasta -t mRNA -o HAN412_Eugene_curated_v1_1_mRNA.fa
#agat_sp_extract_sequences.pl -g HAN412_Eugene_curated_v1_1.gff3 -f Ha412HOv2.0-20181130.fasta -t ncRNA -o HAN412_Eugene_curated_v1_1_ncRNA.fa
#cat HAN412_Eugene_curated_v1_1_mRNA.fa HAN412_Eugene_curated_v1_1_rRNA.fa HAN412_Eugene_curated_v1_1_ncRNA.fa > HAN412_Eugene_curated_v1_1_all_transcripts.fa

# extract just the longest isoform for each Trinity 'gene'
python2 ~/gsd_RNA-seq/scripts/diff_iso/get_longest_isoforms.py \
	~/gsd_RNA-seq/data2/transcriptome/all_sample_Trinity.cd-hit-est_99.single_line.fasta \
	> ~/gsd_RNA-seq/data2/transcriptome/longest_isos_list.all_sample_Trinity.cd-hit-est_99.txt

cut -f1 longest_iso_list.all_sample_Trinity.cd-hit-est_99.txt \
	| grep -A1 --no-group-separator -Ff - ~/gsd_RNA-seq/data2/transcriptome/all_sample_Trinity.cd-hit-est_99.single_line.fasta \
	> ~/gsd_RNA-seq/data2/transcriptome/all_sample_Trinity.cd-hit-est_99.longest_isos.single_line.fasta

# BLASTN de novo transcripts against Ha412HO genes
#makeblastdb -in HAN412_Eugene_curated_v1_1_all_transcripts.fa -dbtype nucl
blastn \
	-query ~/gsd_RNA-seq/data2/transcriptome/all_sample_Trinity.cd-hit-est_99.longest_isos.single_line.fasta \
	-db ~/gsd_RNA-seq/data/ref_genome_Ha412HO/HAN412_Eugene_curated_v1_1_all_transcripts.fa \
	-outfmt 6 \
	> ~/gsd_RNA-seq/analysis/BLAST_out/Trinity_transcripts_vs_HAN412_transcripts.blast

# also run blast to output single top hit
blastn \
	-query ~/gsd_RNA-seq/data2/transcriptome/all_sample_Trinity.cd-hit-est_99.longest_isos.single_line.fasta \
	-db ~/gsd_RNA-seq/data/ref_genome_Ha412HO/HAN412_Eugene_curated_v1_1_all_transcripts.fa \
	-max_target_seqs 1 \
	-outfmt 6 \
	> ~/gsd_RNA-seq/analysis/BLAST_out/Trinity_transcripts_vs_HAN412_transcripts.top_hits.blast

cut -f1,2 ~/gsd_RNA-seq/analysis/BLAST_out/Trinity_transcripts_vs_HAN412_transcripts.top_hits.blast | sort -u > temp
cut -f1 temp | sed 's/_i.*//g' > temp2
paste temp2 temp > ~/gsd_RNA-seq/analysis/BLAST_out/Trinity_transcripts_vs_HAN412_transcripts.top_hits.txt


