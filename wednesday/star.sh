#!/bin/bash

#SBATCH --array=1-24  # NEED TO CHANGE THIS!
#SBATCH --job-name=star # Job name
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=360
#SBATCH --mem=7000 # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --output=arrayJob_%A_%a.out # File to which STDOUT will be written
#SBATCH --error=arrayJob_%A_%a.err # File to which STDERR will be written
#SBATCH --reservation=workshop


begin=`date +%s`

echo $HOSTNAME
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

sample=`sed "${SLURM_ARRAY_TASK_ID}q;d" samples.txt`
R1=${sample}.sickle.R1.fastq
R2=${sample}.sickle.R2.fastq

mkdir ${sample}_star_alignment
REFDIR=
GTF=

module load star

STAR --runThreadN 4 \
--sjdbOverhang 99 \
--genomeDir $REFDIR \
--sjdbGTFtagExonParentGene gene_id \
--sjdbGTFfile $GTF \
--outSAMtype BAM Unsorted SortedByCoordinate \
--outReadsUnmapped Fastx \
--quantMode GeneCounts \
--outFileNamePrefix ${sample}_star_alignment/${sample}_ \
--readFilesIn $R1 $R2 \
> ${sample}_star_alignment/${sample}_STAR.stdout 2> ${sample}_star_alignment/${sample}_STAR.stderr

end=`date +%s`
elapsed=`expr $end - $begin`
echo Time taken: $elapsed

