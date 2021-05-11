#!/bin/zsh -l

#SBATCH --time=14-00:00:00   # Use the form DD-HH:MM:SS
#SBATCH --mem=40G # Memory Per node
#SBATCH --mail-user=lianggao@ksu.edu --mail-type=ALL   # same as =BEGIN,FAIL,END
#SBATCH --job-name=bcftools.Aet.owwc.275.258.937
#SBATCH --nodes=1
#SBATCH --cpus-per-task=25
#SBATCH --partition=ksu-plantpath-jpoland.q,batch.q,ksu-gen-highmem.q
#SBATCH --array=1-7
#SBATCH --output="log/log_%x_%a_%A.txt"

cores="$SLURM_CPUS_PER_TASK"; jobname="$SLURM_JOB_NAME"
jobid="$SLURM_ARRAY_JOB_ID"; taskid="$SLURM_ARRAY_TASK_ID"


prefix="200626.Aetv4.owwc.bcftools.${taskid}"
ref='/bulk/jpoland/genome/Ae_tauschii/AL8_78/index/GCA_002575655.1_Aet_v4.0_genomic.fna'
crams="/homes/lianggao/owwc/data/200626.bam.list.wgs275.gbs258.nextera937.txt" #can be list of CRAM/BAM files
chr_intervals_folder='/homes/lianggao/owwc/data'
intervals=$(ls $chr_intervals_folder/190114.*CM*.txt | sed -n ${taskid}p)

mkdir $prefix


#run SNP calling for the first hundred bins
cat  $intervals | parallel --will-cite -j $cores  /homes/lianggao/scripts/Variant_call/call_bcftools.zsh $ref $crams $prefix '{}' 



