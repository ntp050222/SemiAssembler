#/bip7_disk/Yu_Ting_102/samtools-0.1.19/misc/wgsim -e 0 -d 450 -s 50 -N 2321000 -1 101 -2 101 -r 0 -R 0 -X 0 sequence.fa read_pe1.fq read_pe2.fq
#/bip7_disk/Yu_Ting_102/samtools-0.1.19/misc/wgsim -e 0 -d 4500 -s 50 -N 2321000 -1 101 -2 101 -r 0 -R 0 -X 0 sequence.fa read_mt1.fq read_mt2.fq


ref=ref.fa
read_pe1=read_pe1.fq
read_pe2=read_pe2.fq

read_mt1=read_mt1.fq
read_mt2=read_mt2.fq

bwa index -a bwtsw $ref

bwa aln -n 0 -o 0 -t 20 $ref $read_pe1 > pe_1.sai
bwa aln -n 0 -o 0 -t 20 $ref $read_pe2 > pe_2.sai
bwa sampe $ref pe_1.sai pe_2.sai $read_pe1 $read_pe2 > pe.sam

bwa aln -t 20 $ref $read_mt1 > mt_1.sai
bwa aln -t 20 $ref $read_mt2 > mt_2.sai
bwa sampe $ref mt_1.sai mt_2.sai $read_mt1 $read_mt2 > mt.sam

samtools view -Sb pe.sam > pe.bam
samtools sort pe.bam pe_sort
samtools view -h -o pe_sort.sam pe_sort.bam
samtools index pe_sort.bam

samtools view -Sb mt.sam > mt.bam
samtools sort mt.bam mt_sort
samtools view -h -o mt_sort.sam mt_sort.bam
samtools index mt_sort.bam

rm pe.sam
rm pe_sort.sam
rm mt.sam
rm mt_sort.sam
