
date >> total_time

ref=reference.fa
mv $ref ref_0.fa
read_pe1=Read_pe1.fq
read_pe2=Read_pe2.fq

read_mt1=Read_mt1.fq
read_mt2=Read_mt2.fq

java -jar /bip7_disk/Yu_Ting_102/picard-tools-1.114/NormalizeFasta.jar I=ref_0.fa O=ref.fa


mv $read_pe1 read_pe1.fq
mv $read_pe2 read_pe2.fq
mv $read_mt1 read_mt1.fq
mv $read_mt2 read_mt2.fq



echo "[[Create SAM]]"
./create_samfile.sh
echo "[[Detect insertion and deletion]]"
./detect_ID.py

date >> total_time

echo "[[Insertion assemble]]"
mkdir assemble
cd assemble
ln -s ../sh_sga.sh .
ln -s ../sub.sh .
ln -s ../blast_insertion.sh
ln -s ../find_primer.py .
echo "[[Create Graph]]"
./sh_sga.sh

date >> total_time

echo "[[Assemble sequence]]"
./sub.sh
./blast_insertion.sh
cd ..

date >> total_time

cat ./assemble/output_Assemble_2 output_deletion_range > output_I_D_tn
sort -k1,1 -k2,2n output_I_D_tn > output_I_D_tn_sort
echo "[[Draft genome]]"
./insert_SV.py

./create_contigSAM.sh
./contig_replacement.py

date >> total_time
