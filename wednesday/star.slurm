#!/bin/bash

#SBATCH --job-name=star # Job name
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=60
#SBATCH --mem=20000 # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --partition=gc
#SBATCH --reservation=workshop
#SBATCH --array=1-24
#SBATCH --output=slurmout/arrayJob_%A_%a.out # File to which STDOUT will be written
#SBATCH --error=slurmout/arrayJob_%A_%a.err # File to which STDERR will be written

start=`date +%s`
echo $HOSTNAME
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

sample=`sed "${SLURM_ARRAY_TASK_ID}q;d" samples.txt`
REF=
GTF=

outpath='02-STAR_alignment'
[[ -d ${outpath} ]] || mkdir ${outpath}
[[ -d ${outpath}/${sample} ]] || mkdir ${outpath}/${sample}

echo "SAMPLE: ${sample}"

module load star
call="STAR --runThreadN 8 \
--sjdbOverhang 99 \
--genomeDir $REF \
--sjdbGTFtagExonParentGene gene_id \
--sjdbGTFfile $GTF \
--outSAMtype BAM SortedByCoordinate \
--outReadsUnmapped Fastx \
--quantMode GeneCounts \
--outFileNamePrefix ${outpath}/${sample}/${sample}_ \
--readFilesCommand zcat \
--readFilesIn 01-HTS_Preproc/${sample}/${sample}_R1.fastq.gz 01-HTS_Preproc/${sample}/${sample}_R2.fastq.gz"

echo $call
eval $call

end=`date +%s`
runtime=$((end-start))
echo $runtime
