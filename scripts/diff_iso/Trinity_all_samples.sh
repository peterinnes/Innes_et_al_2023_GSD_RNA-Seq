# run Trinity for de novo transcriptome assembly, using reads from one plant from each population (three dune, three non-dune) as input.
# Input samples are specified in Trinity_samples_file.txt
# Trinity v2.13.2

fqdir=~/gsd_RNA-seq/data/fastq/trimmed/
samples_file=/home/peter/gsd_RNA-seq/Trinity_all_samples_file.txt
outdir=~/gsd_RNA-seq/data2/transcriptome/

cd $fqdir

Trinity --seqType fq --max_memory 100G --samples_file $samples_file --CPU 12 --output ${outdir}all_sample_Trinity --verbose --full_cleanup | tee ${outdir}all_sample_Trinity_stdout.txt

