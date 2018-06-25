#!/bin/bash

## assumes star

start=`date +%s`
echo $HOSTNAME

outpath='02-STAR_alignment'
[[ -d ${outpath} ]] || mkdir ${outpath}

REF=""
GTF=""

or sample in `cat samples.txt`
do
  [[ -d ${outpath}/${sample} ]] || mkdir ${outpath}/${sample}
  echo "SAMPLE: ${sample}"

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
done

end=`date +%s`
runtime=$((end-start))
echo $runtime
