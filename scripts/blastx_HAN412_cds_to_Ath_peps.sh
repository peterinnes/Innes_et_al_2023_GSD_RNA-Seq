blastx -query ~/gsd_RNA-seq/data/ref_genome_Ha412HO/HAN412_Eugene_curated_v1_1_cds.fa -db ~/gsd_RNA-seq/data/ref_genome_Araport11/Araport11_pep_20220103.fa -evalue 1e-20 -num_threads 10 -max_target_seqs 1 -outfmt 6 -out ~/gsd_RNA-seq/analysis/BLAST_out/HAN412_vs_Araport11_1e-20.blastx

cd ~/gsd_RNA-seq/analysis/BLAST_out/ 

cut -f1,2 HAN412_vs_Araport11_1e-20.blastx | sort -u > HAN412_vs_Araport11_1e-20.txt
