

date
date
date


SGAA=sga
THREAD=25
OVERLAPLEM=51
MAXSIZE=101

FASTA1=../read_pe1.fq
FASTA2=../read_pe2.fq
PREFIX=target

echo "[[[[[sga preprocess]]]]]"
$SGAA preprocess -o sga_pre_${PREFIX}.fastq --pe-orphans sga_pre_${PREFIX}_orphans.fastq -p 1 $FASTA1 $FASTA2

echo "[[[[[sga index]]]]]"
$SGAA index -t $THREAD --algorithm ropebwt sga_pre_${PREFIX}.fastq #read length <200bp

echo "[[[[[sga correct]]]]]"
$SGAA correct -t $THREAD --algorithm=overlap -o sga_corr_${PREFIX}.fastq sga_pre_${PREFIX}.fastq

echo "[[[[[sga index]]]]]"
$SGAA index -t $THREAD --algorithm ropebwt sga_corr_${PREFIX}.fastq

echo "[[[[[sga rmdup]]]]]"
$SGAA rmdup -t $THREAD -o sga_rmdup_${PREFIX}.fastq sga_corr_${PREFIX}.fastq


echo "[[[[[sga overlap]]]]]"
$SGAA overlap -t $THREAD -m $OVERLAPLEM --exact sga_rmdup_${PREFIX}.fa

sga assemble -m 55 --transitive-reduction -o contig sga_rmdup_${PREFIX}.asqg.gz
sga index contig-contigs.fa


date
date
date
