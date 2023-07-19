Great Sand Dunes prairie sunflower (*Helianthus petiolaris*) ecotype transcriptomics
====================================

This directory contains code from [Innes et al. (2023)](https://doi.org/10.1101/2023.04.22.537924) for the analysis of RNA-seq data from *H. petiolaris* dune and non-dune ecotypes in Great Sand Dunes National Park, Colorado, USA. Scripts (within scripts/) are listed in analysis order.


##### Seedling traits and sampling map
analyze_seedling_traits.R  
plot_sampling_map.R 
  

##### Trimming and alignment:
fastp_trim.sh  
STAR_genomeGenerate.sh  
STAR_align.sh
  

##### SNP calling:
process_bams_for_gatk.sh  
gatk_haplotype_caller.sh  
filter_gatk_vcf.sh  
count_nHet.sh  
run_SNPeff.sh


##### Fst estimates:
calculate_Fst.sh  
get_gene_windows.sh  
calculate_Fst_per_gene.py


##### Expression quantification:
count_transcripts_htseq-count.sh


##### DESeq2 and rMATS analysis:
analyze_expression_DESeq2.R  
rmats.py #copy of rMATS source code. We changed the filtering scheme on lines 335–339.  
run_rMATS.sh  
analyse_splicing_rMATS.R


##### Smith et al differential splicing pipeline. Adapted from [Smith et al. 2021](https://github.com/chriscrsmith/SunflowerAberrantSplicing). Scripts with same number can be run in either order:
diff_iso/00_Trinity_all_samples.sh  
diff_iso/01_remove_redundant_transcripts_cd-hit.sh  
diff_iso/02_blast_transcriptome_to_genome.sh  
diff_iso/02_blast_Trinity_transcripts_to_HAN412HO_transcripts.sh  
diff_iso/03_analyzeHAblast_v11.py  
diff_iso/04_count_isoforms_RSEM.sh  
diff_iso/05_filter_low_expression.R  
diff_iso/06_verifySplicing_v5.py  
diff_iso/07_make_splicedIso_list.py  
diff_iso/08_consolidate_isoforms_v2.py  
diff_iso/09_post_consolidate_filter.py  
diff_iso/10_parents_diff_v2.py  
diff_iso/11_fdr_correction_parents_diff_results.R  
diff_iso/11_post_processing.sh  
diff_iso/12_analyze_iso_proportions.R  
diff_iso/12_compare_rMATS_and_Smith_et_al.sh  
diff_iso/get_longest_isoforms.py


##### Coexpression analysis:
WGCNA/WGCNA.R  
WGCNA/iterativeWGCNA.R  
WGCNA/analyze_iterativeWGCNA.R


##### GO enrichment analysis:
blastx_HAN412_cds_to_Ath_peps.sh  
parse_Ath_GO_terms.py  
get_Ha412HO_GO_terms.R  
run_GO_enrichment_and_clustering.sh  
analyze_GO_enrichment.R


##### Look for spliceosomal homologs:
blast_spliceosome_to_transcripts.sh


##### Downstream analyses:
pca.R #PCA for figure 2  
analyze_inversions.R #enrichment of DE and DS genes within major inversions  
analyze_overlaps.R #overlaps between gene sets  
analyze_Goebl_2022_loci.R #proximity of DE and DS genes to Goebl et al 2022 adaptive loci  
analyze_Fst.R #Fst of different gene sets, Fst vs LFC/dPSI regressions, etc  


##### CODE FOR MAIN TEXT FIGURES. This supercedes code within other scripts for main text figures:
figures_main_text.R
