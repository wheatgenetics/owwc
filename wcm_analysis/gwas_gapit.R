# load required libraries
# source("http://www.bioconductor.org/biocLite.R")
# BiocManager::install('multtest')
# install.packages("gplots")
# BiocManager::install('snpStats')
# install.packages("LDheatmap")
# install.packages("genetics")
# install.packages("ape")
# install.packages("EMMREML")
# install.packages("scatterplot3d") 
library(data.table)
library(multtest)
library(gplots)
library(LDheatmap)
library(genetics)
library(ape)
library(EMMREML)
library(compiler) #this library is already installed in R
library("scatterplot3d")
# source("http://zzlab.net/GAPIT/gapit_functions.txt")
source('~/scripts/R/functions/GWAS/gapit_functions.txt'); #keep a local copy
# source("http://zzlab.net/GAPIT/emma.txt")
source('~/scripts/R/functions/GWAS/emma.txt'); #keep a local copy

# read pheno file
myY = read.delim('data/wcm_pheno.txt')

# read hap file
myG = read.delim('/bulk/lianggao/owwc/vcf.cmb.l1.l2/190625.all.chr.ipkmm.vcf.format.head.hmp.txt', head = F)

myKI = read.csv('data/GAPIT.Kin.VanRaden.csv', head = F) # kinship matrix
myCV = read.csv('data/GAPIT.PCA.csv',head = T) # principal components

# run gwas
myGAPIT = GAPIT(Y = myY, G = myG, KI = myKI, CV = myCV)

