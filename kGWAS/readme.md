# kGWAS

## Description
kGWAS is an extension of [AgRenSeq_GLM](https://github.com/kgaurav1208/AgRenSeq_GLM), a pipeline to identify candidate resistance (_R_) genes in plants directly from a diversity panel. Please refer to the original pipeline for all the steps. Only the modifications made to the original pipeline are described below. These modifications pertain to the following steps of the original pipeline:

### 1: Split reference fasta file
The script `ref_splitter.py` splits each chromosome of the reference line into a number of equal chunks.  The following command will split each chromosome of al878 into chunks of size 10000:
 
````
python ref_splitter.py -a al878.fa -s 10000 -o al878_split.fa 
````

### 2: Generate association scores of _k_-mers and project onto the split assembly

```
python RunAssociation_GLM.py -i presenceMatrix.txt.gz -hd presenceMatrix_header.txt -a al878_split.fa -p phenotype.txt -s snp.tsv  -u usable.txt -o output.txt
```



## Pre-requisites

### Python 3 and above

The code has been tested in Python 3.5.3. 

The following Python modules are required:

* `numpy` 
* `pandas` 
* `Biopython`: to parse assembly file
* `scikit-learn`: to compute PCA from SNP markers matrix
* `statsmodels`: for regression analysis
* `bitarray`



#### Parameters

Parameter (long version)| Parameter (short version) | Argument | Description
--- | --- | --- | ---
--inputmatrix | -i | presenceMatrix.txt.gz | Mandatory. The path to file containing the gzipped version of presence/absence matrix of k-mers.
--header | -hd | presenceMatrix_header.txt | Mandatory. The path to file containing the list of accessions in the order that their presence is scored in the presence/absence matrix.
--assembly | -a | assembly.fasta | Mandatory. The path to split reference assembly onto which k-mers are mapped for plotting.
--phenotype | -p | phenotype.txt | Mandatory. The path to phenotype file.
--usable | -u | usable.txt | Optional. The path to file containing the list of usable accessions.
--stackman | -st | snp.txt | Mandatory. The path to file containing matrix of SNP markers to compute PCA and correct for population structure.
--subsample | -sub | integer | Optional. Run association analysis by taking a random subsample of _this size_ from the given accessions.
--permute | -per |  | Optional. Permute the phenotype scores.
--pcadimensions | -dim | integer | Default 3. The Number of significant PCA dimensions used as covariates for regression analysis.
--correlationthreshold  | -c | float | Default 0.2. Only those k-mers whose correlation with phenotype is greater than _this value_ are retained for  regression analysis
--pvalthreshold  | -pv | float | Default 15.0. Only those k-mers with pvalue greater than _this value_ are retained
--output | -o | output.txt | Mandatory. The path to output file to store the association values corresponding to nlr contigs in the given assembly.
