# run Trinity for de novo transcriptome assembly of individual samples.
# Trinity v2.13.2

fqdir=~/gsd_RNA-seq/data/fastq/trimmed/
samples_file=/home/peter/gsd_RNA-seq/Trinity_samples_file.txt
outdir=~/gsd_RNA-seq/data/transcriptome/

samples=$(cut -f3 $samples_file | cut -d '_' -f1-2)

cd $fqdir
for sample in $samples; do
	echo "assembling $sample..."
	Trinity --seqType fq --max_memory 50G --single ${sample}_R1_001.trimmed.fq.gz \
		--CPU 8 --output $outdir${sample}_Trinity/ --verbose --full_cleanup \
		| tee ${outdir}${sample}_Trinity_stdout.txt
done
