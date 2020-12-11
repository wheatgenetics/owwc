#!/bin/bash -l

#SBATCH --job-name=vcf2ipkmm.filter.shuff.2long
#SBATCH --partition=ksu-plantpath-jpoland.q,batch.q,killable.q
#SBATCH --mem-per-cpu=4G             # Memory per core/cpu
#SBATCH --cpus-per-task=1
#SBATCH --array=1-1010             # 1-4 or 1,2,3-4 or stepwise 1-5:2 for 1,3,5
#SBATCH --time=00-23:00:00           # DD-HH:MM:SS
#SBATCH --mail-type=ALL
#SBATCH --output="log/log_%x_%a_%A.txt"

## set environmental variables
cores="$SLURM_CPUS_PER_TASK"; jobname="$SLURM_JOB_NAME"
jobid="$SLURM_ARRAY_JOB_ID"; taskid="$SLURM_ARRAY_TASK_ID"

## input files vcf.gz raw
vcf_gz_list='/homes/lianggao/owwc/201026.owwc.wgs.275.ex31.vcf.gz.list.txt'

## path to required software
filter_vcf='/homes/lianggao/scripts/vcf_filtering/filter_vcf_rm_path.awk'
vcf2long_v1='/homes/lianggao/owwc/scripts/vcf2long.v1.sh'

## create required directories
currDir="$SLURM_SUBMIT_DIR"
cd "${currDir}"

vcf_gz=$(cat "${vcf_gz_list}"| sed -n "${taskid}"p)
echo "The vcf name is ${vcf_gz}. Job started on $(date)"; echo
	
filtered_vcf=${vcf_gz/vcf.gz/filtered.vcf}
zcat $vcf_gz  | $filter_vcf -v dphom=2 -v dphet=4 -v minqual=40 -v mindp=100 -v minhomn=1 -v minhomp=0.9 -v tol=0.2 -v minmaf=0.01 -v minpresent=0.8 > $filtered_vcf
echo "VCF filtering  minqual=40 -v mindp=100 -v minhomn=1 -v minhomp=0.9 -v tol=0.2 -v minmaf=0.01 -v minpresent=0.8$(date)"; echo

grep "#" $filtered_vcf > $filtered_vcf.100.shuf
grep -v "#" $filtered_vcf | shuf -n 100 | sort -V -k1 -k2n >> $filtered_vcf.100.shuf
echo "VCF shuffle (random sel of 100?) variants done on $(date)"; echo

cat $filtered_vcf.100.shuf | $vcf2long_v1 | grep -Fv './.'| bgzip -@$cores  > $filtered_vcf.100.shuf.2long.gz
echo "Converted shuffled vcf to long format to do all vs all allele comparisons"


tsv2pct_ident='/homes/lianggao/owwc/scripts/Rscripts/06_geno_pct_ident_vcflong2yy_generic_all.v.all.R'
tsv_gz=$filtered_vcf.100.shuf.2long.gz
percent_ident=${tsv_gz/.gz/.percent.ident.txt}
module load R/3.5.1-foss-2018b ## note need a lib from this, but the actual R version is not 3.5.1
echo "VCF long format convert to percent identity for $tsv_gz started on $(date)"; echo
/homes/lianggao/.local/bin/Rscript  $tsv2pct_ident $tsv_gz $percent_ident
echo "VCF long format convert to percent identity for $tsv_gz ended on $(date)"; echo

