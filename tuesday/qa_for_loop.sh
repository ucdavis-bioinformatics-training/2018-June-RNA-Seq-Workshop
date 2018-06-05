#!/bin/bash

#SBATCH --nodes=1
#SBATCH --time=600
#SBATCH --mem=2000 # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --reservation=workshop


begin=`date +%s`
echo $HOSTNAME

for sample in `cat samples.txt`
do

R1=${sample}_L006_R1_001.fastq.gz
R2=${sample}_L006_R2_001.fastq.gz

module load scythe
module load sickle

scythe -a adap.fasta -q sanger -o ${sample}.scythe.R1.fastq $R1
scythe -a adap.fasta -q sanger -o ${sample}.scythe.R2.fastq $R2
sickle pe -f ${sample}.scythe.R1.fastq -r ${sample}.scythe.R2.fastq -t sanger -o ${sample}.sickle.R1.fastq -p ${sample}.sickle.R2.fastq -s ${sample}.singles.fastq

done

end=`date +%s`
elapsed=`expr $end - $begin`
echo Time taken: $elapsed

