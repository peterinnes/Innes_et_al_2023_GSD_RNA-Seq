# run Trinity for de novo transcriptome assembly, using reads from one plant from each population (three dune, three non-dune) as input.
# Input samples are specified in Trinity_samples_file.txt
# Trinity v2.13.2

fqdir=~/gsd_RNA-seq/data/fastq/trimmed/
samples_file=/home/peter/gsd_RNA-seq/Trinity_samples_file.txt
outdir=~/gsd_RNA-seq/data/transcriptome/

cd $fqdir

Trinity --seqType fq --max_memory 50G --samples_file $samples_file --CPU 10 --output ${outdir}multi-sample_Trinity --verbose --full_cleanup | tee ${outdir}multi-sample_Trinity_stdout.txt

