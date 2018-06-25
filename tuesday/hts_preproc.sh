#!/bin/bash

## assumes htstream is available on the Pathway

start=`date +%s`
echo $HOSTNAME

outpath='01-HTS_Preproc'
[[ -d ${outpath} ]] || mkdir ${outpath}

for sample in `cat samples.txt`
do
  [[ -d ${outpath}/${sample} ]] || mkdir ${outpath}/${sample}
  echo "SAMPLE: ${sample}"

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
done

end=`date +%s`
runtime=$((end-start))
echo $runtime
