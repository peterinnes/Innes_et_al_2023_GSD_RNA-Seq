# Smith et al aka 'parents diff' DS genes with significant blast hits to Ha412HO
cut -f1 significant_parents_diff_iso_results.2023-07-03.txt | sed 's/\..//g' | grep -Ff - ~/gsd_RNA-seq/analysis/BLAST_out/Trinity_transcripts_vs_HAN412_transcripts.blast | awk '$12>100' | cut -f2 | sort -u | sed 's/.*RNA/gene/g' > HAN412HOv2_all_hits_bit100.signigicant_parents_diff_iso_results.2023-07-03.txt

# overlap b/w Smith et al DS genes and rMATS DS genes
grep -Ff ~/gsd_RNA-seq/analysis/GO_analysis/study_DS_rMATS_genes.txt ~/gsd_RNA-seq/analysis/diff_iso_out/2023-07-03/HAN412HOv2_all_hits_bit100.signigicant_parents_diff_iso_results.2023-07-03.txt| wc -l
# = 290

# Smith et al AS genes with significant blast hits to Ha412HO
cut -f1 good_isoforms.post_consolidate_filtered.2023-07-03.txt | sed 's/\..//g' | grep -Ff - ~/gsd_RNA-seq/analysis/BLAST_out/Trinity_transcripts_vs_HAN412_transcripts.blast | awk '$12>100' | cut -f2 | sort -u | sed 's/.*RNA/gene/g' > HAN412HOv2_all_hits_bit100.AS_genes.2023-07-03.txt

# overlap b/w Smith et al AS genes and rMATS AS genes
grep -Ff ~/gsd_RNA-seq/analysis/rMATS/results_2022-07-14/rMATS_AS_genes.txt HAN412HOv2_all_hits_bit100.AS_genes.2023-07-03.txt | wc -l
# = 4539

# Reciprocal best blast hits between Smith et al DS genes and Ha412HO genes
cut -f1 significant_parents_diff_iso_results.2023-07-03.txt | sed 's/\..//g' | grep -Ff - ~/gsd_RNA-seq/analysis/BLAST_out/Trinity_transcripts_vs_HAN412_transcripts.reciprocal_best_hits.txt | sed 's/.*RNA/gene/g' > HAN412HOv2_reciprocal_best_hits.significant_parents_diff_iso_results.2023-07-3.txt

# overlap b/w Smith et al and rMATS considering just reciprocal best hits
cut -f1 ../2023-02-01/HAN412HOv2_reciprocal_best_hits.significant_parents_diff_iso_results.2023-02-05.txt| grep -Ff - ~/gsd_RNA-seq/analysis/GO_analysis/study_DS_rMATS_genes.txt | wc -l
# = 53
