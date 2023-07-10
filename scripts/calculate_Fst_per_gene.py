repl_python()
import allel
import pandas as pd
import numpy as np

print(allel.__version__)

df_samples = pd.read_table("/home/peter/gsd_RNA-seq/vcf_sample_list.txt")
pop1 = 'dune'
pop2 = 'non-dune'

subpops = {
    pop1: df_samples[df_samples.habitat == pop1].index,
    pop2: df_samples[df_samples.habitat == pop2].index,
}

# sample indices within vcf for pop 1
pop1_idx = subpops[pop1]
# sample indices within vcf for pop 2
pop2_idx = subpops[pop2]

# read-in custom gene windows
all_windows = pd.read_table("analysis/Fst/single_gene_windows_5kb_buffer.bed",
  dtype={'win_start':int})
  
chroms = ['Ha412HOChr01', 'Ha412HOChr02', 'Ha412HOChr03', 'Ha412HOChr04',
'Ha412HOChr05', 'Ha412HOChr06', 'Ha412HOChr07', 'Ha412HOChr08', 'Ha412HOChr09',
'Ha412HOChr10', 'Ha412HOChr11', 'Ha412HOChr12', 'Ha412HOChr13', 'Ha412HOChr14',
'Ha412HOChr15', 'Ha412HOChr16', 'Ha412HOChr17']
  
# make empty list
res_a = []
res_b = []
res_c = []
res_d = []

# loop through chromosomes
for chrom in chroms:
  callset = allel.read_vcf(input=
    "/home/peter/gsd_RNA-seq/data/hc_out/dune_non-dune.nHet_filtered.hc.vcf.gz",
    region=chrom)
    
  #sorted(callset.keys())
  
  #callset['variants/CHROM']
  #callset['samples']
  
  gt = allel.GenotypeArray(callset['calldata/GT'])
  
  snps_pos = allel.SortedIndex(callset['variants/POS'])
    
  windows = all_windows.loc[all_windows['chrom'] ==
    chrom].iloc[:,1:].to_numpy()
  
  f, w, n = allel.windowed_weir_cockerham_fst(pos=snps_pos, g=gt,
    subpops=[pop1_idx, pop2_idx], windows=windows)
  
  chrom_arr = np.full(f.shape, chrom) #make an array of current chrom ID 
  
  res_a.append(f) #Fst
  res_b.append(w) #windows
  res_c.append(n) #n variants
  res_d.append(chrom_arr)

# concatenate the lists of arrays into single arrays
res_a = np.concatenate(res_a)
res_b = np.concatenate(res_b)
win_start, win_end = zip(*res_b) #split 2-D array into multiple 1-D arrays
res_c = np.concatenate(res_c)
res_d = np.concatenate(res_d)


results = pd.DataFrame({
  'chrom': res_d,
  'win_start': win_start,
  'win_end': win_end,
  'Fst': res_a,
  'n_vars': res_c})

results.to_csv("analysis/Fst/Fst_single_gene_windows_5kb_buffer.allel.csv", na_rep='NA', index=False)
