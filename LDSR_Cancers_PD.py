## Linkage disequilibrium score regression (LDSR) between cancers and PD

#We have avlaible full summary statistics from  four studies (breast cancer, pancreatic cancer, endometrial cancer and melanoma)
#Summary statistics for LDSR were pre-prepared and containing columns: 

    data_selected <- select("rsID","Allele1","Alelle2", "Zscore", "P_value", N)
    #Z-score was calculate as Effect (divided) divided to Standart error
    data1$Zscore <- data1$Effect/data1$StdErr

#LDSR was performed using  python package (Bulik-Sullivan B et al., 2015)
munge_sumstats.py --out PD.sumstats --sumstats PD_noUKBB.txt --merge-alleles ../w_hm3.snplist --N-col N --snp SNP --a1 A1 --a2 A2 --p Pvalue --signed-sumstats Zscore,0
munge_sumstats.py --out endometrial --sumstats endometrial_cancerSE_N.txt --merge-alleles w_hm3.snplist --N-col N --snp snpid --a1 A1 --a2 A2 --p Pvalue --signed-sumstats Zscore,0
munge_sumstats.py --out brest --sumstats endometrial_cancerSE_N.txt --merge-alleles w_hm3.snplist --N-col N --snp snpid --a1 A1 --a2 A2 --p Pvalue --signed-sumstats Zscore,0
munge_sumstats.py --out Prostate_cancer --sumstats endometrial_cancerSE_N.txt --merge-alleles w_hm3.snplist --N-col N --snp snpid --a1 A1 --a2 A2 --p Pvalue --signed-sumstats Zscore,0
munge_sumstats.py --out melanoma --sumstats Melanoma_GWAS.txt --merge-alleles ../w_hm3.snplist --N-col N --snp SNP --a1 A1 --a2 A2 --p Pvalue --signed-sumstats Zscore,0
munge_sumstats.py --out kerat --sumstats ../keratinocytes/keratinocytes_summstats.txt --merge-alleles ../w_hm3.snplist --N-col N --snp MarkerName --a1 Allele1 --a2 Allele2 --p P-value --signed-sumstats Zscore,0

#LDSR between cancers and PD summary statistics
ldsc.py --rg PD.sumstats.sumstats.gz,Prostate_cancer.sumstats.gz --out PD_Prostate --ref-ld-chr ../eur_w_ld_chr/ --w-ld-chr ../eur_w_ld_chr/
ldsc.py --rg PD.sumstats.sumstats.gz,breast.sumstats.gz --out PD_Breast --ref-ld-chr ../eur_w_ld_chr/ --w-ld-chr ../eur_w_ld_chr/
ldsc.py --rg PD.sumstats.sumstats.gz,endometrial.sumstats.gz. --out PD_Endometrial --ref-ld-chr ../eur_w_ld_chr/ --w-ld-chr ../eur_w_ld_chr
ldsc.py --rg PD.sumstats.gz,melanoma.sumstats.gz --out PD_Melanoma --ref-ld-chr ../eur_w_ld_chr/ --w-ld-chr ../eur_w_ld_chr/
ldsc.py --rg PD.sumstats.sumstats.gz,kerat.sumstats.gz --out PD_kerat --ref-ld-chr ../eur_w_ld_chr/ --w-ld-chr ../eur_w_ld_chr/
