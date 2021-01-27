## Pre-requisites

### Jellyfish
Install [Jellyfish](https://github.com/gmarcais/Jellyfish) with Python binding.

### Python 3 and above

The following Python modules are required:

* `Biopython`: to parse assembly file

## Description

- Count k-mers from either trimmed or raw read files for each Ae. tauschii and wheat landrace accession using jellyfish.

```
zcat accession1_R?.fastq.gz | jellyfish count -C -m 51 -s 3G -o accession1.jf /dev/fd/0
```

- Create configuration files - `lin1_jf.cfg`, `lin2_jf.cfg`, `lin3_jf.cfg`, `landrace_jf.cfg` -  containing the paths to the jellyfish of each Ae. tauchii lineage 1 accession, Ae. tauchii lineage 2 accession, Ae. tauchii lineage 3 accession and wheat landrace, respectively. For example, `lin1_jf.cfg` would look like:

```
path/to/lineage1_accession1.jf
path/to/lineage1_accession2.jf
...
```

- Split D subgenome of wheat assembly into scaffolds of desired window size (eg: 100000)  using script `ref_splitter.py` from kGWAS:
 
````
python ref_splitter.py -a arina_subgenome_D.fa -s 100000 -o arinaD_split.fa 
````

- Chop the output from the previous step into files of approximately equal size for parallel processing

````
python chop_assembly.py -a arinaD_split.fa -c 100000000 -o arinaD_chop
````

- Run the following script for each assembly chunk in parallel:

```
python admixture.py -a arinaD_chop_0.fa -l1c lin1_jf.cfg -l2c lin2_jf.cfg -l3c lin3_jf.cfg -lrc landrace_jf.cfg -o output_0.txt
```

`-lrc`  is an optional parameter and not required for assemblies of wheat landraces, such as, Chinese Spring.
