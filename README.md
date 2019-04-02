# YeastRVM
## Bloom
This folder has some supplementary materials from Bloom et al.'s paper: <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4001867/>  
 
**nature11867-s2.xls** - It contains drug doses, heritability statistics, and QTL summary statistics for traits investigated in their study.  
**nature11867-s3.xls** - This table shows the additive genetic variance, partitioned by chromosome, for each trait.  
**nature11867-s4.xls** - It containes detected QTL, positions, effect sizes, confidence intervals, and genes underneath detected QTL for each trait.  
**BYxRM_GenoData.csv** - Genotype data for each segregant as a tab-delimited text file.


## Grinberg
This folder has some files from Grinberg et al.'s paper: <https://www.escholar.manchester.ac.uk/uk-ac-man-scw:269442>

**cross.RData** - Original phenotypic and genotypic data of Bloom et al.  
**geno_1.rds** - Unaltered genoset used in their study

## Data
**geno.csv** - It contains the yeast genotypes: 1,008 yeast strains of 11,623 SNPs (extracted from **cross.RData**)  
**pheno.csv** - It contains the yeast phenotypes: <= 46 traits for each strain (extracted from **geno_1.rds**)

## Kernels
**kernel_cosine** - Cosine similarity kernel  
**kernel_rbf1.csv** - Gaussian kernel with *gamma=1e-4*  
**kernel_rbf2.csv** - Gaussian kernel with *gamma=2e-4*  
**kernel_rbf3.csv** - Gaussian kernel with *gamma=5e-5*  
**kernel_rbf4.csv** - Gaussian kernel with *gamma=8.6e-5*  
**kernel_rbf5.csv** - Gaussian kernel with *gamma=3e-4*  **kernel_3gram.csv** - 3-gram kernel  **kernel_5gram.csv** - 5-gram kernel   **kernel_10gram.csv** - 10-gram kernel

## Results
**gsr\_yeast1\_effSNPs.csv** - Ensemble RVM Results: This table shows the number of hits for each SNP for Trait 1 (Cadmium Chloride).   

**gsr\_yeast20\_effSNPs.csv** - Ensemble RVM Results: This table shows the number of hits for each SNP for Trait 20 (Lithium Chloride).   

**gsr\_yeast24\_effSNPs.csv** - Ensemble RVM Results: This table shows the number of hits for each SNP for Trait 24 (Mannose).  

**gsr\_yeast\_CV\_all.csv** - This table contains the cross validation results of single RVM predictions in all traits. 
 
- Columns:  
	1. **trait\_no** : 1-46  
	2. **try\_no** : the m-th try of cross validation (here, 1-10)  
	3. **kernel\_name** : (1:cosine similarity, 2:rbf1, 22:rbf2, 23:rbf3, 24:rbf4, 25:rbf5, 31:3gram, 32:5gram 33:10gram)
	4. **nfolds** : the number of folds in cross validation (here, 10)   
	5. **R2_test** : Calculated R2 (coefficient of determination) for the test fold  

**yeast_result.xlsx** : Summary of cross validation results, Comparisons (with Grinberg's), and heritibility calculations
	

## Scripts
**GSr_RVM.m** - A MATLAB function which uses [SparseBayes][] to implement kernel RVM for predicting quantitative traits 

**GSr_FS.m** - A MATLAB function which uses [SparseBayes][] to implement ensemble of linear basis RVMs for identifying influential markers on a trait 

[SparseBayes]: http://www.miketipping.com/sparsebayes.htm.