Intro to the Command Line
==========================

Joe Fass

jnfass@ucdavis.edu


Download the slides [here](CLIntro.pdf)

Multi-line commands from slides:
--------------------------------

    srun -t 1440 -n 1 --mem 8000 --reservation workshop --pty /bin/bash

    wget --no-check-certificate ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/PhiX/Illumina/RTA/PhiX_Illumina_RTA.tar.gz
    curl -k https://bioshare.bioinformatics.ucdavis.edu/bioshare/download/pdhqicmfgw2bra8/variant.neighborhoods.fa > regions.fa
   
    grep ">" regions.fa | cut -c 2-
    grep ">" regions.fa | cut -c 2- | cut -f1 -d:
    grep ">" regions.fa | cut -c 2- | cut -f1 -d: | sort
    grep ">" regions.fa | cut -c 2- | cut -f1 -d: | sort | uniq -c
    grep ">" regions.fa | cut -c 2- | cut -f1 -d: | sort | uniq -c | sort -rn -k1,1

    ln -s PhiX/Illumina/RTA/Sequence/WholeGenomeFasta/genome.fa .
    
    bwa mem genome.fa regions.fa 1> aln.sam 2> aln.err

    grep ">" regions.fa | cut -c2- | while read header; do echo "contig_$header" >> b; done




 
