# Repository for scripts and data analysis protocols for wheat curl mite (wcm).

Authors: Dr. Liangliang Gao and Paula Silva (KSU, USA)

## Software dependencies 
- SLURM High performance computing job submission system
- HISAT2 V2.1.0
- SAMTOOLS v1.9 
- BCFTOOLS v1.6 
- R 3.6.0 and packages specified in scripts

## 1. Alignment to in-silico reference genome (AB-Jagger + D-tauschii) using hisat2
	aln_wheat_tauschii_cmc.sh (Partial adaptation from Dr. Narinder Singh, KSU, USA)

## 2. Call variants for short arm of chr 6D  using bcftools
	call.snps.cmc.chr6DS.sh
	
## 3. Filter variants using bcftools (step 1) and custom pipeline (step 2)
	Refer to snp_filtering.txt for indications
	vcf_MAF_Missing_Het_filter.pl - script author Dr. Sandesh Shrestha (https://github.com/sandeshsth/01_VCF/tree/master/01_vcf_missing_het_filter, commit log 675f5e1)

## 4. Use the R scripts to perform specific analyses.
	Data folder contains files with required data
	a. gwas_gapit.R - run gwas analysis using GAPIT in R
	b. cluster.R - run hierarchical cluster analysis in R
