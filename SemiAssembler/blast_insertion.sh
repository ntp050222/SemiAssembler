awk '{count++; if(count%5==1 || count%5==2) print $0}'> samread_r1.fa ../LR_read_info
awk '{count++; if(count%5==3 || count%5==4) print $0}'> samread_r2.fa ../LR_read_info

makeblastdb -in Assemble_inserions.fa -dbtype nucl -hash_index -parse_seqids -out Assemble_inserions_db
blastn -query samread_r1.fa -db Assemble_inserions_db -task blastn -dust no -outfmt 6 -out samread_r1.blast
blastn -query samread_r2.fa -db Assemble_inserions_db -task blastn -dust no -outfmt 6 -out samread_r2.blast
#sed 's/_/\t/g' samread_r1.blast > samread_r1.tab
#sed 's/_/\t/g' samread_r2.blast > samread_r2.tab

awk '{if(NR==1){get1=$1;print $0;}else if(NR>1){ if($1==$2 && get1 != $1){print $0;get1=$1}}}' samread_r1.blast > samread_r1.blast.filter
awk '{if(NR==1){get1=$1;print $0;}else if(NR>1){ if($1==$2 && get1 != $1){print $0;get1=$1}}}' samread_r2.blast > samread_r2.blast.filter
#awk '{if($1==$2)print $0}' samread_r2.blast > samread_r2.blast.filter

samtools faidx Assemble_inserions.fa
./find_primer.py

rm samread_*
rm Assemble_inserions_db*
