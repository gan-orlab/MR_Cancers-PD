## Mendelian randomization between cancers and Parkinson's disease
For the construction of our IVs, we selected studies from the GWAS Catalog (Buniello A et al., 2019) using the R package “MRInstruments” which is a part of TwoSampleMR R package (Hemani G et al, 2017; Hemani G et al., 2018). 
Studies were selected based on specific criteria. First, we searched traits by the keywords “cancer”, “carcinoma”, “glioma”, “lymphoma”, “leukemia”, “melanoma”. We then selected the most recent available GWAS for each cancer, with a minimum of 1000 cases and European ancestry. Subsequently, we included SNPs with GWAS significate _p_-values (<5×10−8) from each study.  The final code for the study selection represented below. Additionally, recent GWASs on melanoma and combined analysis of keratinocyte cancers were added as they were not available in the GWAS Catalog.
`#studies selected from the GWAS catalog
devtools::install_github("MRCIEU/MRInstruments")
library(MRInstruments)
data(gwas_catalog)
head(gwas_catalog)
setwd("C:/R/MR/Oncology/Selection/")
#Renal cell carcinoma
renal_carcinoma <- subset(gwas_catalog, grepl("Scelo", Author) & Phenotype == "Renal cell carcinoma" & pval <= 5e-08)
renal_carcinoma_dat <- format_data(renal_carcinoma)
write.csv(renal_carcinoma_dat,file="renal_carcinoma_GWAS.csv",quote = FALSE)
#Breast cancer
breast_carcinoma <- subset(gwas_catalog, grepl("Michailidou", Author) & Phenotype == "Breast cancer" & pval <= 5e-08)
breast_carcinoma_dat <- format_data(breast_carcinoma)
write.csv(breast_carcinoma_dat,file="Breast_cancer_GWAS.csv",quote = FALSE)
#Prostate cancer
prostate <- subset(gwas_catalog, grepl("Schumacher", Author) & Year == 2018 & pval < 5e-08)
prostate_cancer_dat <- format_data(prostate)
write.csv(prostate_cancer_dat,file="prostate_cancer.csv",quote = FALSE)
#Skin cacners
	#cutaneous squamous cell carcinoma (Chahal) - failed to calculate SE
Cutaneous_Squamous_cell_carcinoma <- subset(gwas_catalog, grepl("Chahal", Author) & pval <= 5e-08)
Cutaneous_Squamous_cell_carcinoma_dat <- format_data(Cutaneous_Squamous_cell_carcinoma)
	#SE calculations here and below if absent from GWAS catalog data
Cutaneous_Squamous_cell_carcinoma_dat$se.exposure<- TwoSampleMR:::get_se(Cutaneous_Squamous_cell_carcinoma_dat$beta.exposure,Cutaneous_Squamous_cell_carcinoma_dat$pval.exposure)
write.csv(Cutaneous_Squamous_cell_carcinoma_dat,file="Cutaneous_Squamous_cell_carcinoma.csv",quote = FALSE)
	#Non-melanoma skin cancer
Nonmelanoma_skincancer <- subset(gwas_catalog, grepl("Visconti", Author) & pval <= 5e-08 & Phenotype == "Non-melanoma skin cancer")
Nonmelanoma_skincancer_dat <- format_data(Nonmelanoma_skincancer)
write.csv(Nonmelanoma_skincancer_dat,file="Nonmelanoma_skincancer_dat.csv",quote = FALSE)
	#melanoma - #, abcent from the catalog, IVs created manually from the paper
#Blood cancers
	#Lymphoma (Genome-wide association study of classical Hodgkin lymphoma identifies key regulators of disease susceptibility)
lymphoma_1 <- subset(gwas_catalog, grepl("Sud", Author) & Year == 2017 & pval <= 5e-08)
lymphoma_dat <- format_data(lymphoma_1)
write.csv(lymphoma_dat,file="lymphoma.csv",quote = FALSE)
	#Chronic lymphocytic leukemia
lymphoma_law <- subset(gwas_catalog, grepl("Law P", Author) & Year == 2017 & pval <= 5e-08 & Phenotype == "Chronic lymphocytic leukemia")
lymphoma_law_dat <- format_data(lymphoma_law)
lymphoma_law_dat$se.exposure<- TwoSampleMR:::get_se(lymphoma_law_dat$beta.exposure,lymphoma_law_dat$pval.exposure)
write.csv(lymphoma_law_dat,file="lymphoma_bcell.csv",quote = FALSE)
#Pancreatic cancer
pancreatic <- subset(gwas_catalog, grepl("Klein", Author) &  pval <= 5e-08  & Phenotype == "Pancreatic cancer")
pancreatic_dat <- format_data(pancreatic)
write.csv(pancreatic_dat,file="pancreatic.csv",quote = FALSE)
#Urinary bladder cancer
Urinary_bladder_cancer <- subset(gwas_catalog, grepl("Rafnar", Author) &  pval <= 5e-08 & Year == 2018)
Urinary_bladder_cancer_dat <- format_data(Urinary_bladder_cancer)
write.csv(Urinary_bladder_cancer_dat,file="Urinary_bladder_cancer.csv",quote = FALSE)
#Gastrointestinal system
	#Colorectal cancer
Colorectal_cancer <- subset(gwas_catalog, grepl("Law", Author) &  pval <= 5e-08  & Phenotype == "Colorectal cancer")
Colorectal_cancer_dat <- format_data(Colorectal_cancer)
write.csv(Colorectal_cancer_dat,file="Colorectal_cancer.csv",quote = FALSE)
	#Upper digestive
Upper_digestive <- subset(gwas_catalog, grepl("Lesseur", Author) &  pval <= 5e-08 )
Upper_digestive_dat <- format_data(Upper_digestive)
Upper_digestive_dat$se.exposure<- TwoSampleMR:::get_se(Upper_digestive_dat$beta.exposure,Upper_digestive_dat$pval.exposure)
write.csv(Upper_digestive_dat,file="Upper_digestive.csv",quote = FALSE)
#Glioma/non-glioma 
Glioma <- subset(gwas_catalog, grepl("Melin", Author) &  pval <= 5e-08)
Glioma_dat <- format_data(Glioma)
write.csv(Glioma_dat,file="Glioma_dat.csv",quote = FALSE)
#lung
Lung_cancer <- subset(gwas_catalog, grepl("McKay", Author) & pval <= 5e-08 & Phenotype_simple == "Lung cancer")
Lung_cancer_dat <- format_data(Lung_cancer)
write.csv(Lung_cancer_dat,file="Lung_cancer.csv",quote = FALSE)
#Combined analysis of keratinocyte cancers identifies novel genome-wide loci- table added manually, abcent from the catalog, IVs created manually 
#add  in final script from other page
#Epithelian ovarian cancer
endometrial_cancer <- subset(gwas_catalog, grepl("O'Mara TA", Author) & pval <= 5e-08)
endometrial_cancer_dat <- format_data(endometrial_cancer)
write.csv(endometrial_cancer_dat,file="endometrial_cancer.csv",quote = FALSE)quote = FALSE)
#Cervical cancer 
Cervical_cancer <-subset(gwas_catalog,grepl("Cervical cancer", Phenotype) & pval <= 5e-08 & Year == 2017)
Cervical_cancer_dat <- format_data(Cervical_cancer)
Cervical_cancer_dat$se.exposure<- TwoSampleMR:::get_se(Cervical_cancer_dat$beta.exposur,Cervical_cancer_dat$pval.exposure)
write.csv(Cervical_cancer_dat,file="Cervical_cancer.csv",quote = FALSE)`

## R2 and F-statistics was calculated for each study
R2- proportion of variance in exposure variable explained by SNPs 
F-statistics ‘strength’ of the genetic instrumental variable

    rFromN <-get_r_from_lor(exp_data$beta.exposure, exp_data$eaf.exposure, exp_data$ncase.exposure, xp_data$ncontrol.exposure ,0.01,  model = "logit")
    mean(rFromN)
    R2 <- mean(rFromN)
    R2
    #F-statistics
	#k = number of variants in the study
	#n = sample size
	F<-(R2*(n-1-k))/((1-R2)*k)
	F

## Mendelian randomization between all cancers as exposure and PD as outcome (PD summary statistics without the UK biobank cohort)

We created for loops, to run all MR analyses simultaneously. For that, we created folders for each study.
After we enterted R and activated TwoSampleMR and ggplot2. All 

    #creating folders for each study
    for chr in Breast_cancer Colorectal_cancer Cutaneous_Squamous_cell_carcinoma endometrial_cancer glioma keratinocytes lungcancer lymphoma_bcell lymphoma melanoma Nonmelanoma_skincancer pancreaticcancer prostate_cancer renal_carcinoma Upper_digestive Urinary_bladder_cancer; do mkdir ${chr} ; done
    #Open R, install packages if nessesary.
    install.package("TwoSampleMR")
    require(TwoSampleMR)
    install.package("ggplot2")
    require(ggplot2)
    #add all exposures in the loop
    #this, allow us to run all MR analyses simultaniously, saving to the appropriate folders
    file <- c ("Breast_cancer_clumped.csv","Colorectal_cancer_clumped.csv","Cutaneous_Squamous_cell_carcinoma_clumped.csv", "endometrial_cancer_clumped.csv", "glioma_clumped.csv", "keratinocytes_clumped.csv", "lungcancer_clumped.csv", "lymphoma_bcell_clumped.csv", "lymphoma_clumped.csv", "melanoma_clumped.csv","Nonmelanoma_skincancer_clumped.csv","pancreaticcancer_afterclumped.csv","prostate_cancer_clumped.csv","renal_cancer_clumped.csv", "Upper_digestive_clumped.csv","Urinary_bladder_cancer_clumped.csv")
    folder <- c("Breast_cancer", "Colorectal_cancer", "Cutaneous_Squamous_cell_carcinoma", "endometrial_cancer", "glioma", "keratinocytes", "lungcancer" , "lymphoma_bcell", "lymphoma", "melanoma", "Nonmelanoma_skincancer", "pancreaticcancer", "prostate_cancer", "renal_carcinoma", "Upper_digestive", "Urinary_bladder_cancer")
    for(i in 1:length(file)){
    #here we simultaniously read exposure data and perform clumping with standart parameters (Clumping window 10,000 kb, R2 cutoff 0.001) 
    exp_data <- read_exposure_data(file[i], sep=",", snp_col = "SNP", beta_col = "beta.exposure",  eaf_col = "eaf.exposure", se_col="se.exposure", effect_allele_col = "effect_allele.exposure", other_allele_col= "other_allele.exposure", pval_col = "pval.exposure", samplesize_col = "N.exposure", clump=FALSE)
    out_data <- read_outcome_data(snps = exp_data$SNP,filename = "kostiya_snp_all.txt", sep="\t",  snp_col = "SNP", beta_col = "Effect", se_col = "StdErr", eaf_col = "Freq1", effect_allele_col = "Allele1", other_allele_col = "Allele2", pval_col = "P-value", ncase_col = "ncase", ncontrol_col = "ncontrol")
    dat <- harmonise_data(exposure_dat=exp_data, outcome_dat=out_data, action=2)
    #Using the output from the  `mr`  function this report will generate a report containing tables and graphs summarising the results. A separate report is produced for each exposure - outcome pair that was analysed.
    mr_report(dat, study = folder[i],output_path = folder[i]) 
    #to save forest plot in high resolution, using ggplot.
    res_single <- mr_singlesnp(dat)
    p5 <- mr_forest_plot(res_single)
    p5[[1]]
    ggsave(p5[[1]], file= "plot.jpg", path=folder[i] , width=7, height=12)
    }
  
## Mendelian randomization between all cancers as exposure and PD as outcome (PD summary statistics without the UK biobank cohort)

#here we repeated all analyses using PD summary statistics with the UK biobank cohort

    require(TwoSampleMR)
    require(ggplot2)
    file <- c ("Breast_cancer_clumped.csv","Colorectal_cancer_clumped.csv","Cutaneous_Squamous_cell_carcinoma_clumped.csv", "endometrial_cancer_clumped.csv", "glioma_clumped.csv", "keratinocytes_clumped.csv", "lungcancer_clumped.csv", "lymphoma_bcell_clumped.csv", "lymphoma_clumped.csv", "melanoma_clumped.csv","Nonmelanoma_skincancer_clumped.csv","pancreaticcancer_afterclumped.csv","prostate_cancer_clumped.csv","renal_cancer_clumped.csv", "Upper_digestive_clumped.csv","Urinary_bladder_cancer_clumped.csv")
    folder <- c("Breast_cancer", "Colorectal_cancer", "Cutaneous_Squamous_cell_carcinoma", "endometrial_cancer", "glioma", "keratinocytes", "lungcancer" , "lymphoma_bcell", "lymphoma", "melanoma", "Nonmelanoma_skincancer", "pancreaticcancer", "prostate_cancer", "renal_carcinoma", "Upper_digestive", "Urinary_bladder_cancer")
    for(i in 1:length(file)){
    exp_data <- read_exposure_data(file[i], sep=",", snp_col = "SNP", beta_col = "beta.exposure",  eaf_col = "eaf.exposure", se_col="se.exposure", effect_allele_col = "effect_allele.exposure", other_allele_col= "other_allele.exposure", pval_col = "pval.exposure", samplesize_col = "N.exposure", clump=TRUE)
    out_data <- read_outcome_data(snps = exp_data$SNP,filename = "META_no23_yesUKBB.txt", sep="\t",  snp_col = "SNP", beta_col = "b", se_col = "se", eaf_col = "freq", effect_allele_col = "A1", other_allele_col = "A2", pval_col = "p", ncase_col = "N_cases", ncontrol_col = "N_controls")
    dat <- harmonise_data(exposure_dat=exp_data, outcome_dat=out_data, action=2)
    mr_report(dat, study = folder[i],output_path = folder[i]) 
    res_single <- mr_singlesnp(dat)
    p5 <- mr_forest_plot(res_single)
    p5[[1]]
    ggsave(p5[[1]], file= "plot.jpg", path=folder[i] , width=7, height=12)
    }
    
## Reverse MR; PD as exposure, melanoma as outcome

    exp_data <- read_exposure_data(PD_exposure.csv, sep=",", snp_col = "SNP", beta_col = "beta.exposure",  eaf_col = "eaf.exposure", se_col="se.exposure", effect_allele_col = "effect_allele.exposure", other_allele_col= "other_allele.exposure", pval_col = "pval.exposure", samplesize_col = "N.exposure", clump=TRUE)
    out_data <- read_outcome_data(snps = exp_data$SNP,filename = "melanoma.final.meta.gwas.txt", sep="\t",  snp_col = "rsID", beta_col = "Estimate_effect",  eaf_col = "EAF", se_col="SE", effect_allele_col = "Allele1", other_allele_col= "Allele2", pval_col = "P_value", ncase_col = "Case_size", ncontrol_col = "Control_size")
    dat <- harmonise_data(exposure_dat=exp_data, outcome_dat=out_data, action=2)
    mr_report(dat, study = folder[i],output_path = folder[i]) 
    res_single <- mr_singlesnp(dat)
    p5 <- mr_forest_plot(res_single)
    p5[[1]]
    ggsave(p5[[1]], file= "melanoma.jpg", width=7, height=12)
