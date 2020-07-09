
# Aligning whole genome deep sequencing (10-30x) and skim sequencing (GBS or Nextera) tauschii samples to reference AL8, call SNPs and calculate percent identical calls for samples potentially same origin.
Authors: Dr. Liangliang Gao (with scripts adapted from Dr. Narinder Singh (KSU, USA)  and Dr. Martin Mascher (IPK Germany))

## Software dependencies
- SLURM High performance computing job submission system
- HISAT2 V2.1.0
- SAMTOOLS v1.9 
- BCFTOOLS v1.9
- GNU Parallel 20160822
- AWK 4.0.2
- R 3.6.0
- Rscript 3.6.0
- R packages (data.table v1.12.8); (openxlsx 4.1.01); (zoo 1.8-6) 

## (1)  alignment using hisat2
alignment_to_ref

## (2) call SNPs by breaking reference into 4Mb intervals and use GNU parallel, bcftools
call.snps.Aet_v4.cmb.owwc.275.gbs258.nextera937.zsh


## (3) convert VCF file to long format and calculate percent identical calls for samples presumably same origin
vcf2long.2pct.ident.275.258.937.sh


