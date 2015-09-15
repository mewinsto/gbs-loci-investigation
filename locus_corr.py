#!/usr/bin/python

import sys
import itertools
from Bio import SeqIO
import numpy as np

infile = sys.argv[1]
outfile = sys.argv[2]
snp_file = sys.argv[3]
ref_file = sys.argv[4]
C_OUT = sys.argv[5]
D_OUT = sys.argv[6]
MAF_OUT = sys.argv[7]
MAF_MAX = sys.argv[8]
final_out = sys.argv[9]
REF = open(ref_file).readlines()
data = open(infile).readlines()
outfile = open(outfile,'w')
snp_file = open(snp_file,'w')

head = data[0].strip().split(" ")
ntaxa = int(head[1])
loci = int(head[3])
snps = int(head[5])

#DEFINE REF ALLELES
ref = REF[0].lstrip().rstrip().split(" ")
ref = filter(None, ref)
ref = [x for x in ref if x != "_"]

#MAKE SNP-LOCUS INDEX
#index_1 = data[1].lstrip().rstrip().split(" ")
#index_1 = filter(None, index_1)
#index_1 = index_1[1:]
#index_length = len(index_1)
#icount = 1
#snp_index = []
#for loc in index_1:
#    it = [i for i in itertools.repeat(icount,len(loc))]
#    snp_index.append(it)
#    icount += 1
#snp_index = list(itertools.chain.from_iterable(snp_index))
#for item in snp_index:
#    print >>snp_file, item    
    
#MAKE ALLELES BINARY
for line in data[1:]:
    a = line.lstrip().rstrip().split(" ")
    a = filter(None, a)
    a = [x for x in a if x != "_"]
    sample = a[0]
    bin = []
    for x in range(1,len(a)-1):
        locus = a[x]
        z = 0
        while z < len(locus):
            if locus[z] == "N":
                bin.append("-")
            elif locus[z] == ref[x-1][z]:
                bin.append("0")
            else:
                bin.append("1")   
            z += 1    
        bin.append(" ")
    bin = ''.join(bin)

    print >>outfile, bin    
            
outfile.close()

##READ IN OUTFILE AND CALCULATE MAF
##LoL is a list of lists of lists, where the first list is the samples, the second are the loci, and the third are the snps

outfile = sys.argv[2]
binary = open(outfile).readlines() 
LoL = []
for k in range(0,len(binary)):
    e = []
    d = binary[k].lstrip().rstrip().split(" ")
    for i in range(0,len(d)):
        e.append(list(d[i]))
    LoL.append(e)


C_OUT = open(C_OUT,'w')
D_OUT = open(D_OUT,'w')
MAF_OUT = open(MAF_OUT,'w')
MAF_MAX = open(MAF_MAX,'w')
final_out = open(final_out,'w')
d_length = len(d)
LoL_length = len(LoL)
het_max_c = np.empty([LoL_length,1])          
for j in range(0,d_length):
        tmp = []
        locus_length = len(LoL[0][j])
        for i in range(0,LoL_length):
                tmp.append(LoL[i][j])
        Y = np.array(tmp).reshape(LoL_length,locus_length)
        c_list = []
        d_list = []
        maf_snp = []
        for i in range(0,locus_length):
                baba = Y[:,i]
                C = 0
                D = 0
                for k in range(0,len(baba)):              
                    if baba[k] == "0":
                        C += 1
                    elif baba[k] == "1":
                        D += 1
                c_list.append(C)
                d_list.append(D)        
                C = float(C)
                D = float(D)
                maf = D/(C+D)
                maf_snp.append(maf)

        c_list_s = [str(x) for x in c_list]
        d_list_s = [str(x) for x in d_list]
        maf_snp_s = [str(x) for x in maf_snp]
        c_out = ' '.join(c_list_s)
        d_out = ' '.join(d_list_s)
        maf_out = ' '.join(maf_snp_s)
        C_OUT.write(c_out)
        C_OUT.write("\n")
        D_OUT.write(d_out)
        D_OUT.write("\n")
        MAF_OUT.write(maf_out)
        MAF_OUT.write("\n")
        MAF_MAX.write(str(max(maf_snp)))
        MAF_MAX.write("\n")
        
        for i in range(0,len(c_list)):
            if maf_snp[i] == max(maf_snp):
                het_max = Y[:,i]
        het_max_c = np.column_stack((het_max_c,het_max))

het_max_final = het_max_c[:,1:]
for i in range(0,LoL_length):
        samp_bin = ' '.join(het_max_final[i,:])
        final_out.write(samp_bin)
        final_out.write("\n")
#        print >>big_test_output, samp_bin
         
C_OUT.close()
D_OUT.close()
MAF_OUT.close()
MAF_MAX.close()
final_out.close()
