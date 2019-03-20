---
title: 'Paternal BMI PACE Analysis Plan: Phase One'
author: "Gemma Sharp"
date: "02/01/2018"
output:
  word_document:
    fig_caption: yes
    toc: yes
    toc_depth: 2
  pdf_document:
    fig_caption: yes
    toc: yes
    toc_depth: 2
  html_document:
    fig_caption: yes
    toc: yes
    toc_depth: 2
---

# Background and aims
There is increasing evidence that paternal lifestyle factors and environmental exposures can influence offspring health, potentially via epigenetic changes to the germline. Previous studies have provided some evidence that paternal obesity is associated with offspring DNA methylation at imprinted genes involved in early growth and metabolism. Other potential explanations for associations between paternal obesity and offspring methylation include unobserved and/or residual confounding, including that introduced by maternal adiposity.

The aims of this PACE project are to investigate:

1. Associations between paternal periconceptional obesity and offspring blood DNA methylation at birth, childhood and beyond.

2. The likelihood that such associations are causal, i.e. due to a direct effect of paternal obesity rather than a maternal (intrauterine) effect or confounding by genetic variation or other pre- or post-natal environmental factors.

3. Offspring phenotypes that might be mediated by paternal obesity-associated DNA methylation.

The project is split into two phases. **The first phase, described in this analysis plan, will focus on aim 1**. 

Findings from the first phase will feed into the second phase, where we will contact cohorts again to complete additional analyses relating to aims 2 and 3. In aim 2, we will use three causal inference techniques, the maternal/paternal negative control design, cross-cultural comparisons and Mendelian randomization, in an attempt to infer causality in associations identified in aim 1. In aim 3, we will use mQTLs (SNPs robustly associated with methylation) and MR-Base (a database and analytical platform for Mendelian randomization) to conduct a pheWAS (phenome-wide association study) to identify offspring phenotypes associated with differential methylation at paternal obesity-associated CpGs.

All of the code necessary to run these analyses is provided here https://github.com/ammegandchips/PACE_Paternal_BMI/blob/master/Rcode.R. We hope that this will lessen the burden for contributing cohorts and it will also ensure that all outputs are consistently formatted. This will significantly reduce the burden on the people conducting the meta-analysis. Therefore, **please use the supplied code**.

# Analyses relating to aim 1

## Exposure

* **pat.bmi**: Paternal BMI in kg/m2, either self-reported, maternal-reported or measured. Please double check BMI values >= +/- 5 SD from the mean in your dataset to make sure they are not data entry errors.

* **mat.bmi**: In order to compare paternal and maternal effects, some models include maternal BMI (as either the main exposure or a covariate). Pre-pregnancy BMI is preferred, but early pregnancy BMI is also acceptable. Values can be measured or self-reported and should be in kg/m2. Please double check BMI values >= +/- 5 SD from the mean in your dataset to make sure they are not data entry errors.

* For the EWAS, paternal and maternal BMI will be converted to Z-scores, but *please note that the code will do this for you*.

## Outcome

* Illumina Infinium 450k or EPIC BeadChip DNA methylation data in blood.

* We are interested in newborns (cord blood or neonatal blood spots collected at birth), young children (6 months-4 years old), children (5-11 years old) adolescents (12-17 years old) and/or adults (18 and above). If you have data for multiple time points, please analyse these separately. 

* The methylation data should be normalised beta values on a scale of 0 to 1 with no transformation (e.g. not M values). You can use your preferred method to normalise the data, but our preference is for Functional Normalisation. Please contact [gemma.sharp@bristol.ac.uk](gemma.sharp@bristol.ac.uk) if you would like R code to conduct functional normalisation on your data.

* Outliers should be trimmed using the IQR3 (Tukey) method. The code for doing this is provided.

* Please use your preferred study QC settings for probe filtering. However, please do not exclude probes just because they are on a published list of possibly problematic probes (e.g. Chen or Naeem) and please do not exclude probes on the sex chromosomes. If in any doubt, please include rather than exclude probes at this stage.

## Other variables

Perhaps moreso than in previous PACE analyses, it is very important that covariates are coded exactly as outlined below. The R code relies on these codings!

* **Paternal social class (ses):** binary numeric variable (1=high/0=low), please use your preferred classification, but note that our preference is for education level. *If your study does not have information on paternal social class, or there are a lot of missing values, please use maternal social class instead and note this in your Excel output.*

* **Paternal age (pat.age):** continuous numeric variable in years

* **Maternal age (mat.age):** continuous numeric variable in years

* **Paternal smoking status around pregnancy (pat.active.smoking):** binary numeric variable (1=smoking during pregnancy or <=3 months before conception/0=no smoking during pregnancy or <= 3 months before conception) 

* **Maternal smoking status during pregnancy (mat.active.smoking):** binary numeric variable (1=smoking throughout pregnancy/0=no smoking in pregnancy or quitting smoking after the first trimester)

* **Parity (parity):** binary numeric variable (1=one or more previous children/0=no previous children)

* **Surrogate variables (SVs) to adjust for batch:** please do NOT include a known batch variable in your models or adjust for batch using another method such as ComBat. The code for calculating surrogate variables is encorporated in the EWAS code provided. We hope (with some support for this hope from the literature and personal experience) that using this approach in all cohorts will reduce heterogeneity and lambdas.

* **Estimated cell proportions:** Cell proportions are estimated using the Houseman method (e.g. by using the estimateCellCounts() function in minfi). Studies with newborn methylation should use the Bakulski et al. cord blood reference panel and *include all 7 cell types* generated: nRBC,CD8T, CD4T, NK, Bcell, Mono, Gran. Studies of older children and adults should use the Reinius adult reference panel and *include all 6 cell types* generated: CD8T, CD4T, NK, Bcell, Mono, Gran. 

* **Selection factors:** Please include if relevant for your study, for example if your sample contains cases and controls for a condition, please include the case/control variable (coded as appropriate). If you want to adjust for a selection factor, you will have to add the column name(s) to the objects called traits.and.covariates and covariates, set at lines 77 and 78 of the analysis code.

* **Ethnicity:** If your study has more than one major ethnic group (for example, European ancestry, Latino, African Ancestry, Asian), please analyse them separately.

* **Child's sex (sex):** Binary numeric variable. This will be used to stratify analyses (1=females,0=males).

## Exclusions

* Please exclude data for partners who are not the biological father of the child, according to either maternal or partner report. If there is any uncertainty about whether the partner is the biological father, please exclude. Please make a note of how many fathers are excluded for this reason and include this information in the Excel output.

* Please also exclude multiple pregnancies (e.g. twins) and siblings (i.e. each mother/father should appear only once in the dataset)

## EWAS models

### Minimally-adjusted: Paternal BMI + Cells (model name: min.pat)
Methylation ~ Paternal BMI + SVs for batch + Estimated cell counts

### Minimally-adjusted: Maternal BMI + Cells (model name: min.mat)
Methylation ~ Maternal BMI + SVs for batch + Estimated cell counts

### Minimally-adjusted: adjusted for other parent's BMI and cells (model name: min.mutual)

*Note that the code extracts results for the effect of paternal AND maternal BMI*

Methylation ~ Paternal BMI + Maternal BMI + SVs for batch  + Estimated cell counts

### Covariates-adjusted: Paternal BMI + Covariates + Cells (model name: covs.pat)
Methylation ~ Paternal BMI + SVs for batch + Paternal age + Paternal smoking status + Paternal socioeconomic status + Maternal age + Maternal smoking status + Parity + Estimated cell counts

### Covariates-adjusted: Maternal BMI + Covariates + Cells (model name: covs.mat)
Methylation ~ Maternal BMI + SVs for batch + Paternal age + Paternal smoking status + Paternal socioeconomic status + Maternal age + Maternal smoking status + Parity + Estimated cell counts

### Covariates-adjusted: adjusted for other parent's BMI, covariates and cells (model name: covs.mutual)

*Note that the code extracts results for the effect of paternal AND maternal BMI*

Methylation ~ Paternal BMI + Maternal BMI + SVs for batch + Paternal age + Paternal smoking status + Paternal socioeconomic status + Maternal age + Maternal smoking status + Parity + Estimated cell counts

### Mutually-adjusted, in boys only (model name: mutual.boys)

*Note that the code extracts results for the effect of paternal AND maternal BMI*

*The R code will remove female offspring*

Methylation ~ Paternal BMI + Maternal BMI + SVs for batch + Paternal age + Paternal smoking status + Paternal socioeconomic status + Maternal age + Maternal smoking status + Parity + Estimated cell counts

### Mutually-adjusted, in girls only (model name: mutual.girls)

*Note that the code extracts results for the effect of paternal AND maternal BMI*

*The R code will remove male offspring*

Methylation ~ Paternal BMI + Maternal BMI + SVs for batch + Paternal age + Paternal smoking status + Paternal socioeconomic status + Maternal age + Maternal smoking status + Parity + Estimated cell counts

# Outputs

Please supply the following files in the specified formats:

1) EWAS results (For each time point, one Rdata file labelled as YOURSTUDY.patbmi.ewasresults.timepoint.Rdata. This will contain all outputs from all EWAS for one time point. If you have multiple time points, you will have an Rdata file for each time point: birth, early_childhood, late_childhood, adolescence, or adult)
2) IQR log file: for each time point, one Rdata file labelled as YOURSTUDY.patbmi.logIQR.timepoint.Rdata.The code will generate this file for you. If you have multiple time points, you will have an Rdata file for each time point: birth, early.childhood, late.childhood, adolescence, or adult.
3) Extra cohort information: one Excel file labelled as YOURSTUDY.patbmi.cohortinfo.xlsx. This will contain information for each cohort relating to things like normalisation method and number of fathers excluded because of non-paternity. Please download and use the template available at [https://github.com/ammegandchips/PACE_Paternal_BMI/blob/master/YOURSTUDY.patbmi.cohortinfo.xlsx](https://github.com/ammegandchips/PACE_Paternal_BMI/blob/master/YOURSTUDY.patbmi.cohortinfo.xlsx)
4) EWAS variables summary: for the main EWAS models (min.mutual and covs.mutual), a csv file labelled as YOURSTUDY.patbmi.modelname.summary.timepoint.csv. This will contain summary statistics summarising predictor variables. The code will generate these files for you. If you have multiple time points, you will have two csv files for each time point: birth, early.childhood, late.childhood, adolescence, or adult.
5) Cell type summary statistics: a csv file labelled as YOURSTUDY.patbmi.cells.res.summary.timepoint.csv. This will contain summary statistics for linear regressions of cell types ~ paternal BMI. If you have multiple time points, you will have one csv files for each time point: birth, early.childhood, late.childhood, adolescence, or adult.

# Upload, timescale and contacts

* The deadline for upload of results is 1st March 2018.
* When you are ready to upload your results, please email [gemma.sharp@bristol.ac.uk](gemma.sharp@bristol.ac.uk) and I will provide you with a personal URL for upload.

# R code

All the R code to perform these analyses is provided at:  [https://github.com/ammegandchips/PACE_Paternal_BMI/blob/master/Rcode.R](https://github.com/ammegandchips/PACE_Paternal_BMI/blob/master/Rcode.R)

**Please** use this code! If you have any questions about the analysis and/or are struggling to get the code to run, please email [gemma.sharp@bristol.ac.uk](gemma.sharp@bristol.ac.uk).

If you have insufficient data to complete one or more of the EWAS, you can just skip those models.
The code also produces .csv files summarising the variables included in the EWASs.
You shouldn't have to rewrite or add to the code, unless otherwise stated.

There are just two inputs required for these analyses:

1) pheno: a dataframe containing all the "phenotype" data needed for this project. Each row is a sample(individual) and each column is a different variable. Necessary variable names are: "pat.bmi", "mat.bmi", pat.active.smoking", "mat.active.smoking", "sex","ses", "pat.age", "mat.age", "parity". If these columns are named differently in your dataset, please rename the columns accordingly.

2) meth: a matrix of methylation illumina beta values. Each column is a sample and each row is a probe on the array (450k or EPIC). Column names must correspond to the sample.id column in pheno.

# Thank you!

Finally, thank you VERY much for taking the time to run these analyses. I will keep you updated on the progress of the project through the PACE calls. If we find associations in Phase One, I will be in touch regarding additional analyses for Phase Two :)
