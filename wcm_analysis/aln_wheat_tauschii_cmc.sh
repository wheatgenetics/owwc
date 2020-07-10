#!/bin/bash -l

#SBATCH --job-name=wheat_tauschii_cmc
#SBATCH --partition=ksu-plantpath-jpoland.q,batch.q,ksu-gen-highmem.q
#SBATCH --gres=killable:0
#SBATCH --mem-per-cpu=6G             # Memory per core/cpu
#SBATCH --cpus-per-task=6
#SBATCH --array=1-92               # 1-4 or 1,2,3-4 or stepwise 1-5:2 for 1,3,5
#SBATCH --time=03-00:00:00           # DD-HH:MM:SS
#SBATCH --mail-type=ALL
#SBATCH --output="log/log_%x_%a_%A.txt"

## ============================================================================
## bash scripts to align wheat_tauschii with or without cmc gene sample raw reads to wheat or wheat_tauschii hybrid reference, sort bam, rmdup...

## path to required software
hisat2='/homes/lianggao/software/hisat2-2.1.0/hisat2'
samtools='/homes/jpoland/tools/samtools/bin/samtools'


## set environmental variables
cores="$SLURM_CPUS_PER_TASK"; jobname="$SLURM_JOB_NAME"
jobid="$SLURM_ARRAY_JOB_ID"; taskid="$SLURM_ARRAY_TASK_ID"

## create required directories
currDir="$SLURM_SUBMIT_DIR"
cd "${currDir}"
mkdir -p bam bam_temp

## Ae. tauschii ref genome
	dbPath="/bulk/jpoland/genome/wheat-tauschii/index/jagger.v1.1.rm.D.tauschii.AL8.chr1to7.fa"

## WGS data
	seq_files="/bulk/mpsilva/wcm/aln_wheat_tauschii_insilicosyn/aln_200513/data_raw2"

## get sample name (same as the directory name containing WGS)
## assign sample name to the array task id
	sample=$(ls "${seq_files}" | sed -n "${taskid}"p)
	echo "The sample name is ${sample}. Job started on $(date)"; echo
		
## convert the list of fastq files to comma separated values and assign to variables
## forward reads
	forward=$(ls "${seq_files}"/"${sample}"/*R1.fq* | tr '\n' ',' | sed 's/,*$//g')
## reverse reads
	reverse=$(ls "${seq_files}"/"${sample}"/*R2*.fq* | tr '\n' ',' | sed 's/,*$//g')

## aligning reads to ref genome using bowtie2 mem
	echo; echo "=== ALIGNMENT and convert to BAM STARTED AT $(date) ==="; echo
	$hisat2 -p "${cores}" -x "${dbPath}" -1 "${forward}" -2 "${reverse}" | $samtools view -@ "${cores}" -f 2 -bhS - > bam/"${sample}".bam 
	echo; echo "=== ALIGNMENT DONE AND BAM OUTPUT COMPLETED AT $(date) ==="; echo
## samtools sort and remove duplicates
## sort and remove duplicates
	$samtools sort -@ "${cores}" -T bam_temp/"${sample}"  bam/"${sample}".bam | $samtools rmdup -S - bam/"${sample}"_sort_rmdup.bam
	echo; echo "=== BAM SORT DONE AND RM DUPLICATES COMPLETED AT $(date) ==="; echo
## index bam file
    $samtools index -c bam/"${sample}"_sort_rmdup.bam

## rename standard output files
	mv log/log_"${jobname}"_"${taskid}"_"${jobid}".txt log/log_"${jobname}"_"${taskid}"_"${jobid}"_"${sample}".txt
	echo; echo "=== LOG FILE RENAMING DONE ==="; echo

## print out a success messagamtools sort -@ "${cores}" -T bam_temp/"${sample}"  bam/"${sample}".bam | $samtools rmdup -S - bam/"${sample}"_rmdup.bam
	echo "========================================================="
	echo "=== Computation finished at $(date) with error code $? ==="
	echo "========================================================="