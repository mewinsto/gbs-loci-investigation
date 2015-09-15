##READ IN OUTFILE AND CALCULATE MAF
##LoL is a list of lists of lists, where the first list is the samples, the second are the loci, and the third are the snps
import numpy as np

big_test = open('SNPs_20cov_out').readlines() 
LoL = []
for k in range(0,len(big_test)):
    e = []
    d = big_test[k].lstrip().rstrip().split(" ")
    for i in range(0,len(d)):
        e.append(list(d[i]))
    LoL.append(e)


C_OUT_20cov = open("C_OUT_20cov",'w')
D_OUT_20cov = open("D_OUT_20cov",'w')
MAF_OUT_20cov = open("MAF_OUT_20cov",'w')
MAF_MAX_20cov = open("MAF_MAX_20cov",'w')
big_test_20cov_output = open("big_test_20cov_output",'w')
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
        C_OUT_20cov.write(c_out)
        C_OUT_20cov.write("\n")
        D_OUT_20cov.write(d_out)
        D_OUT_20cov.write("\n")
        MAF_OUT_20cov.write(maf_out)
        MAF_OUT_20cov.write("\n")
        MAF_MAX_20cov.write(str(max(maf_snp)))
        MAF_MAX_20cov.write("\n")
#        print >>C_OUT, c_out
#        print >>D_OUT, d_out
#        print >>MAF_OUT, maf_out
        
        for i in range(0,len(c_list)):
            if maf_snp[i] == max(maf_snp):
                het_max = Y[:,i]
        het_max_c = np.column_stack((het_max_c,het_max))

het_max_final = het_max_c[:,1:]
for i in range(0,LoL_length):
        samp_bin = ' '.join(het_max_final[i,:])
        big_test_20cov_output.write(samp_bin)
        big_test_20cov_output.write("\n")
#        print >>big_test_output, samp_bin
         
C_OUT_20cov.close()
D_OUT_20cov.close()
MAF_OUT_20cov.close()
MAF_MAX_20cov.close()
big_test_20cov_output.close()
