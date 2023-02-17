transcriptome=~/gsd_RNA-seq/data2/transcriptome/all_sample_Trinity.cd-hit-est_99
fqdir=~/gsd_RNA-seq/data/fastq/trimmed/
outdir=~/gsd_RNA-seq/data2/RSEM_out_23-2-1/

#extract-transcript-to-gene-map-from-trinity $transcriptome.fasta $transcriptome.iso_gene_map.txt

#rsem-prepare-reference \
#    --bowtie2 \
#    --transcript-to-gene-map $transcriptome.iso_gene_map.txt \
#    $transcriptome.fasta \
#    $transcriptome

cd $fqdir
ls *.fq.gz | while read fq; do
    sample_number=$(echo $fq | cut -d '_' -f1 | cut -d '-' -f3)
    if [[ $sample_number -le 15 ]];
        then sample_name=dune_$sample_number;
        else sample_name=non-dune_$sample_number;
    fi;
    
    echo running RSEM on $fq, sample is $sample_name...;

    rsem-calculate-expression \
        $fq \
        $transcriptome \
        $outdir$sample_name \
        --num-threads 20 \
		--strandedness reverse \
		--bowtie2
done

rsem-generate-data-matrix `ls $outdir*.isoform*results* | sort -V | paste -sd' '` > ${outdir}isoform_counts_matrix.txt
rsem-generate-data-matrix `ls $outdir*genes.results | sort -V | paste -sd' '` > ${outdir}gene_counts_matrix.txt
