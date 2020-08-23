# MR_Cancers-PD
In the current study, we applied two-sample MR framework to examine whether certain types of cancers have causal relationship with Parkinson's disease (PD).
In this project, we used GWAS significant SNPs that are associated with a certain cancer as instrumental variables (IVs) in a second GWAS of PD.   

1. For the construction of our IVs, we selected studies from the GWAS Catalog using the R package “MRInstruments” which is a part of TwoSampleMR R package.
2. Overall, 16 cancer studies summary statistics were selected for analysis as exposures. 
3. We calculated F-statistics and R2 for each of the studies in exposure.
4. For the outcome, we used PD GWAS full summary statisctics with and without UK Biobank cohort.
5. We performed for loops, simultaneously with all exposured and one outcome.
6. In each loop, we performed clumping for exposures, harmonization of datasets.
7. Using the output from the 'mr' function, we generated a report containing tables and graphs summarising the results. A separate report is produced for each exposure - outcome pair that was analysed.
8. Report included an inverse-variance weighted (IVW) method, MR-Egger method, and Weighted median (WM).
9. Heterogeneity was tested using Cochran’s Q test in the IVW and MR-Egger methods.
10. For each method, we constructed funnel plots to be able to manually inspect results and notice pleiotropic outliers.
11. We repeated all analyses with PD summary statistics including and excluding UK Biobank cohort to avoid overlap in the samples.
12. Additionally, we performed reverse MR using PD as exposure and Melanoma full summary statistics as outcome.
13. Furthermore, we performed Linlinkage Disequilibrium Score Regression (LDSR) with the four studies (breast cancer, pancreatic cancer, endometrial cancer and melanoma) we have avlaible full summary statistics from. Which was repeated including and excluding UK biobank cohort from PD summary statistics.

The MR_Cancer_PD.R is a script containing all the code we utilized in our study. 
All scripts related to MR are intended for Rstudio. Script related to LDSR is intended for python. Table formating was performed either in R or in Linux.
