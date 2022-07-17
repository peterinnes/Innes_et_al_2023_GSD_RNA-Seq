# process sorted bam files for GATK HaplotypeCaller
# based on GATK RNA-seq best practices: https://github.com/gatk-workflows/gatk4-rnaseq-germline-snps-indels

ref=~/gsd_RNA-seq/data/ref_genome_Ha412HO/Ha412HOv2.0-20181130.fasta
bam_dir=~/gsd_RNA-seq/data/STAR_bams/

#echo "creating reference .dict file"
#gatk CreateSequenceDictionary \
#	    -R $ref

cd $bam_dir
for bam in *Aligned.sortedByCoord.out.bam.gz; do
	sample=$(echo $bam | cut -d "_" -f1-2)
		
	# add read groups
	echo "adding read groups to $sample"
	gatk AddOrReplaceReadGroups \
        	-I $bam \
        	-O $sample.sorted_RGadded.bam \
        	--RGLB lib1 --RGPL illumina --RGPU unit1 --RGSM $sample
    	echo "done adding read group to $sample"

	echo "marking duplicates in $sample..."
	gatk MarkDuplicates \
        	-I $sample.sorted_RGadded.bam \
        	-O $sample.sorted_RGadded_dupmarked.bam \
        	-M $sample.dupmarked_metrics.txt
	
	echo "splitting N cigar reads in $sample..."
	gatk SplitNCigarReads \
        	-R $ref \
        	-I $sample.sorted_RGadded_dupmarked.bam \
        	-O $sample.sorted_RGadded_dupmarked_split.bam
	
	# clean up intermediate bams
	rm $sample.sorted_RGadded.bam
	rm $sample.sorted_RGadded_dupmarked.bam
done
