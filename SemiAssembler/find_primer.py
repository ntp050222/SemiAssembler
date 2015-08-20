#!/usr/bin/env python2.7.9

import sys;
import os;
import time
import pysam
import numpy as np

def ReverseComplement(seq):
    # too lazy to construct the dictionary manually, use a dict comprehension
    seq1 = 'ATCGTAGCatcgtagc'
    seq_dict = { seq1[i]:seq1[i+4] for i in range(16) if i < 4 or 8<=i<12 }
    return "".join([seq_dict[base] for base in reversed(seq)])

def get_chrom(fasta_fh, chrom):
    return fasta_fh.fetch(chrom)

t_start=time.time();
time_list= open("time_findinfo",'w')

output_info= open("output_info_500", 'w')
output_assemble= open("output_Assemble", 'w')
fastafile = "Assemble_inserions.fa"

output_assemble_2= open("output_Assemble_2", 'w')

f=file("samread_r1.blast.filter",'r')
for line in f:
    items=line.strip().split('\t')
    names=items[0].strip().split('_')
    RA_chr=names[0]
    insertion_site=names[1]
    ref_pos_S=names[2]
    ref_pos_E=names[3]
    RA_name=items[1]
    RA_posA=int(items[8])
    RA_posB=int(items[9])
  #  pe_start_pos=int(pos)-100
    f_B=file("samread_r2.blast.filter",'r')
    for line_B in f_B:
        itemsB=line_B.strip().split('\t')
        RB_name=itemsB[1]
        RB_posA=int(itemsB[8])
        RB_posB=int(itemsB[9])
        if(RA_name == RB_name):
            if(RB_posA < RA_posA):  # reverse = 1, need reverse complement
                flag_reverse=1
                posLen=RA_posA-RB_posA
                flag=1
            else:
                flag_reverse=0
                posLen=RB_posA-RA_posA
                flag=2
            
            if(posLen < 10000):
                if(flag == 1):
                    if( RB_posA < RB_posB):
                        POS_S=RB_posA
                    else:
                        POS_S=RB_posB
                    if( RA_posA < RA_posB):
                        POS_E=RA_posB
                    else:
                        POS_E=RA_posA
                elif(flag == 2):
                    if( RA_posA < RA_posB):
                        POS_S=RA_posA
                    else:
                        POS_S=RA_posB
                    if( RB_posA < RB_posB):
                        POS_E=RB_posB
                    else:
                        POS_E=RB_posA
                fasta = pysam.Fastafile(fastafile)
                r_seq = fasta.fetch(RA_name, POS_S-1, POS_E)
                if(flag_reverse == 1):
                    req2=ReverseComplement(r_seq);
                if(flag_reverse == 0):
                    req2=r_seq;
                output_info.write(RA_chr+'\t'+str(insertion_site)+'\t'+str(POS_S)+'\t'+str(POS_E)+'\t'+str(flag_reverse)+'\n');
                output_assemble.write(RA_chr+'\t'+str(insertion_site)+'\t'+ref_pos_S+'\t'+ref_pos_E+'\n'+req2+'\n');
                output_assemble_2.write(RA_chr+'\t'+ref_pos_S+'\t'+ref_pos_E+'\tI\t'+req2+'\n');
    f_B.close()
f.close()
