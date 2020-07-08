#!/bin/zsh


bcftools='/homes/lianggao/miniconda2/envs/selene/bin/bcftools'
tabix='/homes/lianggao/miniconda2/envs/selene/bin/tabix'
bgzip='/homes/lianggao/miniconda2/envs/selene/bin/bgzip'


ref=$1
shift
bamlist=$1
shift
prefix=$1
shift
part=$1


base=$prefix/${prefix}_${part:t}


{

 # record allelic depth for each genotype call (-a DP,DV) 
 # increase per-file read depth (required for merged BAM/CRAM files with many samples)

 $bcftools mpileup -d 1000000 -q 20 -a DP,DV -b $bamlist -r $part -f $ref \
  | $bcftools call -mv -f GQ - | $bgzip > $base.vcf.gz 

 # tabix -C is required for long pseudomolecules (> 530 Mb)

 echo $pipestatus | tee ${base}_pipestatus.txt  | tr ' ' '\n' | grep -q '^[^0$]' || $tabix -Cp vcf $base.vcf.gz 
} 2> $base.err

