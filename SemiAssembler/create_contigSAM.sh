ref=target_s1
read_contig=./assemble/contig-contigs.fa

bwa index -a bwtsw $ref
bwa mem -t 20 $ref $read_contig > contig.sam


samtools view -Sb contig.sam > contig.bam
samtools sort contig.bam contig_sort
samtools view -h -o contig_sort.sam contig_sort.bam
samtools index contig_sort.bam

rm contig.sam
rm contig_sort.sam
