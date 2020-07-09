#!/bin/zsh -l

#SBATCH --time=00-05:00:00   # Use the form DD-HH:MM:SS
#SBATCH --mem=60G # Memory Per node
#SBATCH --mail-user=mpsilva@ksu.edu --mail-type=ALL   # same as =BEGIN,FAIL,END
#SBATCH --job-name=bcftools.wcm.wheat.tauschii
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=30
#SBATCH --partition=ksu-plantpath-jpoland.q,batch.q,ksu-gen-highmem.q

cd /bulk/mpsilva/wcm/aln_wheat_tauschii_insilicosyn/results
bcftools='/homes/lianggao/miniconda3/bin/bcftools'
ref='/bulk/jpoland/genome/wheat-tauschii/index/jagger.v1.1.rm.D.tauschii.AL8.chr1to7.fa'
bams='200516.wheat.tauschii.insilico.bamlist.txt'
output_vcf='200516.0-220Mb.CM008373.1.vcf'
$bcftools mpileup -q 30 -a DP,DV  -b $bams -r "CM008373.1:1-220000000"  -f $ref | $bcftools call -mv -f GQ - > $output_vcf
