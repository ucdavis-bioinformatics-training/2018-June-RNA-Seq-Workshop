#!/bin/bash

#SBATCH --job-name=htstream # Job name
#SBATCH --nodes=1
#SBATCH --ntasks=9
#SBATCH --time=60
#SBATCH --mem=3000 # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --partition=gc
#SBATCH --reservation=workshop
#SBATCH --array=1-24
#SBATCH --output=slurmout/htstream_%A_%a.out # File to which STDOUT will be written
#SBATCH --error=slurmout/htstream_%A_%a.err # File to which STDERR will be written

start=`date +%s`
echo $HOSTNAME
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

sample=`sed "${SLURM_ARRAY_TASK_ID}q;d" samples.txt`

outpath='01-HTS_Preproc'
[[ -d ${outpath} ]] || mkdir ${outpath}
[[ -d ${outpath}/${sample} ]] || mkdir ${outpath}/${sample}

echo "SAMPLE: ${sample}"

module load htstream

call="hts_Stats -O -L ${outpath}/${sample}/${sample}_htsStats.log -1 00-RawData/${sample}/*R1* -2 00-RawData/${sample}/*R2* | \
      hts_SeqScreener -S -O -A -L ${outpath}/${sample}/${sample}_htsStats.log | \
      hts_SuperDeduper -e 250000 -S -O -A -L ${outpath}/${sample}/${sample}_htsStats.log | \
      hts_SeqScreener -s ref/rrna.fasta -r -S -O -A -L ${outpath}/${sample}/${sample}_htsStats.log | \
      hts_AdapterTrimmer -n -S -O -A -L ${outpath}/${sample}/${sample}_htsStats.log | \
      hts_QWindowTrim -n -S -O -A -L ${outpath}/${sample}/${sample}_htsStats.log | \
      hts_NTrimmer -n -S -O -A -L ${outpath}/${sample}/${sample}_htsStats.log | \
      hts_CutTrim -n -m 50 -S -O -A -L ${outpath}/${sample}/${sample}_htsStats.log | \
      hts_Stats -S -A -L ${outpath}/${sample}/${sample}_htsStats.log -g -p ${outpath}/${sample}/${sample}"

echo $call
eval $call

end=`date +%s`
runtime=$((end-start))
echo $runtime
