#!/usr/bin/env python2.7.9

import sys;
import os;
import time
import pysam
import numpy as np


def get_chrom(fasta_fh, chrom):
    return fasta_fh.fetch(chrom)

ref_FILE="ref.fa"

infile_SV = open("output_I_D_tn_sort", 'r')
target = open("target_s1", 'w')

fasta = pysam.Fastafile(ref_FILE)

flag=1;
i=0;
for line in infile_SV:
    i+=1;
    items=line.strip().split('\t')
    chr=items[0]
    pos_start=int(items[1])
    pos_end=int(items[2])-1
    sv_flag=items[3]
    seq=items[4]
    
    flag_1=pos_start-1
    
    if(i==1):
        chr_pre=chr;
        chr_curr=chr;
        curr_chrom_length= fasta.get_reference_length(chr_curr)
        target.write('>'+chr_curr+'\n')
        r_seq = fasta.fetch(chr_curr, 0, flag_1)
        target.write(r_seq)
        if(sv_flag == "I"):
            flag=pos_end
            target.write(seq)
        if(sv_flag == "D"):
            flag=pos_end
    if(i>1):
        chr_curr=chr;
        if(chr_curr != chr_pre):
            r_seq = fasta.fetch(chr_pre, flag, curr_chrom_length)
            target.write(r_seq)
            curr_chrom_length= fasta.get_reference_length(chr_curr)
            chr_pre=chr_curr;
            target.write("\n")
            target.write('>'+chr_curr+'\n')
            r_seq = fasta.fetch(chr_curr, 0, flag_1)
            target.write(r_seq)
            if(sv_flag == "I"):
                flag=pos_end
                target.write(seq)
            if(sv_flag == "D"):
                flag=pos_end
        elif((chr_curr == chr_pre) and (flag < flag_1)):
            r_seq = fasta.fetch(chr_curr, flag, flag_1)
            target.write(r_seq)
            if(sv_flag == "I"):
                flag=pos_end
                target.write(seq)
            if(sv_flag == "D"):
                flag=pos_end
r_seq = fasta.fetch(chr_curr, flag, curr_chrom_length)
target.write(r_seq)
target.close()

