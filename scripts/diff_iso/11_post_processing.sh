#### parse consolidated isoforms for analysis in R #####
cd ~/gsd_RNA-seq/analysis/diff_iso/
cut -f2 consolidated_isoforms.txt | sort -u | grep -v '|' > tmp
cut -f2 consolidated_isoforms.txt | sort -u | grep '|' > tmp2
sed 's/|/\n/g' tmp2 > tmp3
cat tmp tmp3 > consolidated_isoform_clusters.txt
sed 's/,/_/g' consolidated_isoform_clusters.txt | cut -d '_' -f5,10,15,20,25,30,35 > consolidated_isoform_clusters.rownames.txt
paste consolidated_isoform_clusters.rownames.2023-07-03.txt consolidated_isoform_clusters.txt > consolidated_isoform_clusters.with_rownames.txt
rm tmp; rm tmp2; rm tmp3

#### make TPM matrix ####
cd ~/gsd_RNA-seq/data2/RSEM_out_2023-02-01/

# loop through results files and paste them all together
cut -f1 dune_1.isoforms.results > transcript_ids.txt
for file in `ls *isoforms.results | sort --version-sort`; do cut -f6 $file > tmp_$file; done
paste transcript_ids.txt `ls tmp_* | sort -V` > isoform_tpm_matrix.txt

rm tmp*

# make the header
ls *.log | cut -d '.' -f1 | sort -V | tr "\n" "\t" > tmp
echo transcript_id | paste - tmp > header
tail -n +2 isoform_tpm_matrix.txt > tmp
cat header tmp > isoform_tpm_matrix.txt
