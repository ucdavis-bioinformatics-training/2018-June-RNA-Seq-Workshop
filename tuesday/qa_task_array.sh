#!/bin/bash

#SBATCH --array=1-24  # NEED TO CHANGE THIS!
#SBATCH --job-name=trim # Job name
#SBATCH --nodes=1
#SBATCH --time=600
#SBATCH --mem=2000 # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --output=arrayJob_%A_%a.out # File to which STDOUT will be written
#SBATCH --error=arrayJob_%A_%a.err # File to which STDERR will be written
#SBATCH --reservation=workshop


begin=`date +%s`

echo $HOSTNAME
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

sample=`sed "${SLURM_ARRAY_TASK_ID}q;d" samples.txt`
R1=${sample}_L006_R1_001.fastq.gz
R2=${sample}_L006_R2_001.fastq.gz

module load scythe
module load sickle

scythe -a adapters.fasta -q sanger -o ${sample}.scythe.R1.fastq $R1
scythe -a adapters.fasta -q sanger -o ${sample}.scythe.R2.fastq $R2
sickle pe -f ${sample}.scythe.R1.fastq -r ${sample}.scythe.R2.fastq -t sanger -o ${sample}.sickle.R1.fastq -p ${sample}.sickle.R2.fastq -s ${sample}.singles.fastq

end=`date +%s`
elapsed=`expr $end - $begin`
echo Time taken: $elapsed

