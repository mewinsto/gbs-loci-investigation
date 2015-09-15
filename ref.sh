##THIS BASH SRIPT WILL CONSTRUCT THE REFERENCE ALLELES
#!/bin/bash
for i in {13..8574};
do
    cut -b$i burchellii_full_10000min.snps | sort | uniq -c | sort -rn | grep -v 'N' | head -n1 >> ref_count
done

cut -b9 ref_count > REF_SNPs

tr '\n' 'Z' <REF_SNPs> SNPs

sed 's/Z//g' SNPs > SNP_ref_burchellii_full

##ANOTHER OPTION USING AGCT FUNCTION MAY END UP SHORTENING OUTPUT AND CAUSING DOWNSTREAM EFFECTS:
#cut -b$i eciton_SNPs_no_rep_10000min.snps | sort | uniq -c | sort -rn | grep -E 'A|G|C|T' | head -n1 >> ref_count
