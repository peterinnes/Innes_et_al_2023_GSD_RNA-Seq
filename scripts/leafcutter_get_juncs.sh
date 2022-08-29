cd ~/gsd_RNA-seq/data/STAR_bams
for bamfile in `ls *Aligned*.bam`; do
    echo "Converting $bamfile to $bamfile.junc..."
    #samtools index $bamfile
    regtools junctions extract -a 8 -m 50 -s RF -M 500000 $bamfile -o $bamfile.junc
    echo $bamfile.junc >> juncfiles_list.txt
done
