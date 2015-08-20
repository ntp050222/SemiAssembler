#!/bin/bash
THREAD=30
OVERLAPLEM=41
MAXSIZE=101

SGAA=/bip7_disk/Yu_Ting_102/sga_try/weiche_sga/src/SGA/sga

PREFIX=output

target_prefix=contig-contigs
ASQG=contig-graph.asqg.gz


mp_file1=../LR_read_info

rm innerReads.fa outerReads.fa EndReads.fa SeparateReads.fa
rm -r mp-InnerReads mp-OuterReads mp-EndReads mp-SeparateReads
$SGAA subgraph --fasta -o $PREFIX -L 150 -m 1661 -a 200 -d 0 -t 12050 -F 3000000 -V 3000 $target_prefix $mp_file1 $ASQG > stats


cat *-InnerReads.fa > innerReads.fa
cat *-OuterReads.fa > outerReads.fa
cat *-EndReads.fa   > EndReads.fa
cat *-SeparateReads.fa > SeparateReads.fa

sed 's/ /_/g' innerReads.fa    > innerReads_fix.fa
sed 's/ /_/g' outerReads.fa    > outerReads_fix.fa
sed 's/ /_/g' EndReads.fa      > EndReads_fix.fa
sed 's/ /_/g' SeparateReads.fa > SeparateReads_fix.fa



mkdir mp-InnerReads
mv *-InnerReads.fa mp-InnerReads
mkdir mp-OuterReads
mv *-OuterReads.fa mp-OuterReads
mkdir mp-EndReads
mv *-EndReads.fa mp-EndReads
mkdir mp-SeparateReads
mv *-SeparateReads.fa mp-SeparateReads

./blast_insertion.sh

