fqdir=~/data/sunflower/fastq/

for fq in $(ls ${fqdir}raw/); do
       	sample=$(echo $fq | cut -d '.' -f1); echo trimming $sample...

	java -jar ~/software/Trimmomatic-0.39/trimmomatic-0.39.jar \
	SE -threads 14 -phred33 -trimlog ${fqdir}trim.log \
	${fqdir}raw/$fq ${fqdir}trimmed/$sample.trimmed.fq.gz \
	ILLUMINACLIP:/home/peter/software/Trimmomatic-0.39/adapters/TruSeq3-SE.fa:2:20:5:3:true \
	HEADCROP:10 TRAILING:2 MINLEN:65
done
