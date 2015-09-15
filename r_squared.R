##R Script to read and make r^2 matrix

data = read.table("final_out_20covtest", sep = " ", na.strings = "-")
maf = read.table("MAF_MAX_20covtest")
data_matrix = data.matrix(data)
num_snps = length(data_matrix[1,])

r_squared = array(data = NA, dim = c(num_snps,num_snps))

for (i in 1:num_snps){
    for (j in 1:num_snps){
    	if(maf$V1[i] == 0){r_squared[i,j] = 0}
	
	else{if(maf$V1[j] == 0){r_squared[i,j] = 0}    	

	else{r_squared[i,j] = (cor(data_matrix[,i],data_matrix[,j], use = "na.or.complete")^2)}

}
}
}

write(r_squared, file = "/home/mwinston/LOCI_EXPLORE/r_squared_matrix_20covtest", ncolumns = num_snps, sep = " ")