#!/usr/bin/python

import sys
import itertools
from Bio import SeqIO
import numpy as np

infile = sys.argv[1]
outfile = sys.argv[2]
snp_file = sys.argv[3]
ref_file = sys.argv[4]
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
    for x in range(1,1+loci):
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

