# kGWAS

## Description
kGWAS is an extension of [AgRenSeq_GLM](https://github.com/kgaurav1208/AgRenSeq_GLM), a pipeline to identify candidate resistance (_R_) genes in plants directly from a diversity panel. Please refer to the original pipeline for all the steps. Only the modifications made to the original pipeline are described below. These modifications pertain to the following steps of the original pipeline:

### 1: Matrix construction

- Count k-mers from either trimmed or raw read files for each accession using jellyfish.

```
zcat accession1_R?.fastq.gz | jellyfish count -C -m 51 -s 3G -o accession1.jf /dev/fd/0
jellyfish dump -L 2 -ct accession1.jf > accession1.dump.txt
```

- Create a configuration file, `jellies.cfg`, containing the name of the
accession, tab-separated by the path to its jellyfish

```
accession1	path/to/accession1.jf
accession2	path/to/accession2.jf
...
```

- Run the following script for the jellyfish dump file of each accession in parallel:

```
python create_matrix.py -a accession1 -j accession1.dump.txt -c jellies.cfg -o presenceMatrix_accession1.txt
```
Note that the above step requires all the jellyfish files to be present in the working space, but not all the dump files. The dump files can be moved to the secondary storage (eg, hard disks) after generation, brought into the working space only for processing and deleted after processing. Optionally, the output of above script for each accession can be concatenated, sorted and re-split into similar-sized chunks.


### 2: Assembly for k-mer mapping 

Depending on the assembly used to map k-mers, follow either (a) or (b):

#### (a) Split reference fasta file
The script `ref_splitter.py` splits each chromosome of the reference line into a number of equal chunks. The following command will split each chromosome of al878 into chunks of size 10000:
 
````
python ref_splitter.py -a al878.fa -s 10000 -o al878_split.fa 
````


#### (b) Assign coordinates to scaffolds of a non-reference fasta file by anchoring to a reference genome 
Map a non-reference assembly to the reference genome using minimap2 and retain the longest hit:

````
minimap2 -f 0.01 al878.fa bw1111.fa > bw1111_mmap.paf
````

````
awk -v OFS='\t' -F'\t' '{print $1,$4-$3,$6,$8,$9}' bw1111_mmap.paf | sort -k1,1 -k2,2nr |sort -u -k1,1|sort -k3,3 -k4,4n > bw1111_mapping.txt
````

Rename the scaffolds of non-reference assembly to reflect the mapping coordinates with respect to the reference genome:
 
````
python anchor.py -a bw1111.fa -m bw1111_mapping.txt -o bw1111_anchored.fa
````

The next step can be followed either after (a) or (b).
#### (c) Chop assembly file into files of approximately equal size for parallel processing

````
python chop_assembly.py -a bw1111_anchored.fa -c 100000000 -o bw1111_chop
````

### 3: Generate association scores of _k_-mers and project onto the chopped assembly

Run the following command for each chunk of assembly and matrix:

```
python RunAssociation_GLM.py -i presenceMatrix.txt.gz -hd presenceMatrix_header.txt -a bw1111_chop_0.fa -p phenotype.txt -s snp.tsv  -u usable.txt -o output.txt
```

Put all the resultant outputs in a directory and merge them using the following command:

```
python merge_output.py output_dir merged.txt
```

### 4: Bonferroni correction

To estimate the number of multiple tests for Bonferroni correction, you can use the following script:

```
python knum_bonf.py -i presenceMatrix.txt.gz -hd presenceMatrix_header.txt -u usable.txt -o n_kmers.txt
```

### 5: Plotting
The results can be plotted using the following script:

```
python gwas_plot.py merged.txt chromosome_lengths.txt plot.eps
```


## Pre-requisites

### Jellyfish
Install [Jellyfish](https://github.com/gmarcais/Jellyfish) version 2.2.6 or above with Python binding.

### Python 3 and above

The code has been tested in Python 3.5.3. 

The following Python modules are required for GWAS:

* `numpy` (tested with v1.17.4 and v1.18.5)
* `pandas` (tested with v0.23.0 and v1.0.5)
* `Biopython` (tested with v1.72 and v1.77): to parse assembly file
* `scikit-learn` (tested with v0.19.1 and v0.23.1): to compute PCA from SNP markers matrix
* `statsmodels` (tested with v0.9.0. For the more recent versions of `statsmodels`, you might have to change the import statement in `KmerProjection_GLM.py` from `import statsmodels.formula.api as smf` to `import statsmodels.api as smf`): for regression analysis
* `bitarray` (tested with v0.8.1 and v1.4.0)
* `matplotlib` (tested with v3.3.0): for plotting


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
--mincount  | -mc | int | Default 4. Only those k-mers are retained for regression analysis which are present/absent in more than this number of accessions.
--pvalthreshold  | -pv | float | Default 6.0. Only those k-mers with log10 of pvalue greater than _this value_ are retained
--output | -o | output.txt | Mandatory. The path to output file to store the association values corresponding to nlr contigs in the given assembly.
