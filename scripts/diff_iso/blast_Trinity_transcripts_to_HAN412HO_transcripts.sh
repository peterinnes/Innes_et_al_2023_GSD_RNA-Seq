# find specific matches between Trinity 'genes' (just using longest isoform of each 'gene') and annotated HA412HOV2 genes.
cd ~/gsd_RNA-seq/data/ref_genome_Ha412HO/

# conda activate agat
# extract 'biological definition' of mRNAs i.e. exons and UTRs, no introns
#agat_sp_extract_sequences.pl -g HAN412_Eugene_curated_v1_1.gff3 -f Ha412HOv2.0-20181130.fasta -t exon --merge -o HAN412_Eugene_curated_v1_1_exon_merge.fa

# or, this will extract transcript sequences with introns. Use this bc Trinity isoforms could have intronic sequence in the case of intron retention
#agat_sp_extract_sequences.pl -g HAN412_Eugene_curated_v1_1.gff3 -f Ha412HOv2.0-20181130.fasta -t rRNA -o HAN412_Eugene_curated_v1_1_rRNA.fa
#agat_sp_extract_sequences.pl -g HAN412_Eugene_curated_v1_1.gff3 -f Ha412HOv2.0-20181130.fasta -t mRNA -o HAN412_Eugene_curated_v1_1_mRNA.fa
#agat_sp_extract_sequences.pl -g HAN412_Eugene_curated_v1_1.gff3 -f Ha412HOv2.0-20181130.fasta -t ncRNA -o HAN412_Eugene_curated_v1_1_ncRNA.fa
#cat HAN412_Eugene_curated_v1_1_mRNA.fa HAN412_Eugene_curated_v1_1_rRNA.fa HAN412_Eugene_curated_v1_1_ncRNA.fa > HAN412_Eugene_curated_v1_1_all_transcripts.fa


#makeblastdb -in HAN412_Eugene_curated_v1_1_exon_merge.fa -dbtype nucl
#makeblastdb -in HAN412_Eugene_curated_v1_1_all_transcripts.fa -dbtype nucl

blastn -query ~/gsd_RNA-seq/data2/transcriptome/all_sample_Trinity.cd-hit-est_99.longest_isos.single_line.fasta -db ~/gsd_RNA-seq/data/ref_genome_Ha412HO/HAN412_Eugene_curated_v1_1_all_transcripts.fa -outfmt 6 > ~/gsd_RNA-seq/analysis/BLAST_out/Trinity_transcripts_vs_HAN412_transcripts.blast

