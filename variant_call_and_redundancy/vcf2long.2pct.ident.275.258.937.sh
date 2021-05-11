#!/bin/bash -l

#SBATCH --job-name=vcf2long.2pct.ident.275.258.937
#SBATCH --partition=ksu-plantpath-jpoland.q,batch.q,killable.q
#SBATCH --mem-per-cpu=16G             # Memory per core/cpu
#SBATCH --cpus-per-task=1
#SBATCH --array=1-1010              # 1-4 or 1,2,3-4 or stepwise 1-5:2 for 1,3,5
#SBATCH --time=00-23:00:00           # DD-HH:MM:SS
#SBATCH --mail-type=ALL
#SBATCH --output="log/log_%x_%a_%A.txt"

## ============================================================================
## bash scripts to align wgrc gbs demultiplexed reads to AET-V4 ref, sort rmdup.

## path to required software
vcf2long='/homes/lianggao/owwc/scripts/vcf2long.v2.1gt.3dp.4dv.rm.path.awk'


## set environmental variables
cores="$SLURM_CPUS_PER_TASK"; jobname="$SLURM_JOB_NAME"
jobid="$SLURM_ARRAY_JOB_ID"; taskid="$SLURM_ARRAY_TASK_ID"

## create required directories
currDir="$SLURM_SUBMIT_DIR"
cd "${currDir}"

vcf_gz_list='/homes/lianggao/owwc/data/200703.wgs.gbs.nextera.vcf.gz.list.1010.txt'
vcf_gz=$(cat "${vcf_gz_list}"| sed -n "${taskid}"p)
echo "The vcf name is ${vcf_gz}. Job started on $(date)"; echo
	
zcat $vcf_gz  | $vcf2long | grep -Fv './.' |bgzip -@$cores > ${vcf_gz}.tsv.gz   

echo "vcf2long and bgzip compression done on $(date)"; echo


tsv2pct_ident='/homes/lianggao/owwc/scripts/Rscripts/06_geno_pct_ident_vcflong2yy_generic.R'
tsv_gz=${vcf_gz}.tsv.gz
module load R/3.5.1-foss-2018b
echo "tsv2pct_ident for $tsv_gz started on $(date)"; echo
/homes/lianggao/.local/bin/Rscript  $tsv2pct_ident $tsv_gz ${tsv_gz}.yy.txt
echo "tsv2pct_ident completed on $(date)"; echo

