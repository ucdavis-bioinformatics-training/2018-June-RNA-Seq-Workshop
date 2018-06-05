Advanced Command-Line
=======================

**1\.** Advanced Pipes

    mkdir advanced
    cd advanced
    ls /home | cut -c1 | sort | uniq -c
    ln -s /share/biocore-archive/Leveau_J_UCD/RNASeq_Arabidopsis_2016/00-RawData/C61_S67_L006/C61_S67_L006_R*_001.fastq.gz .
    zcat C61_S67_L006_R1_001.fastq.gz | sed -n '1~4p' | cut -d: -f10 | sort | uniq -c

---

**2\.** Process substitution

    module load sickle scythe
    wget https://ucdavis-bioinformatics-training.github.io/2017-June-RNA-Seq-Workshop/tuesday/adapters.fasta
    zcat C61_S67_L006_R1_001.fastq.gz | head -400000 > r1.fq
    sickle se -f <(scythe -a adapters.fasta r1.fq) -t sanger -o trimmed.fa

    zcat C61_S67_L006_R2_001.fastq.gz | head -400000 > r2.fq
    sickle pe -f <(scythe -a adapters.fasta r1.fq) -r <(scythe -a adapters.fasta r2.fq) -t sanger -o trimmed1.fa -p trimmed2.fa -s single.fa

---

**3\.** The 'sed' command

    cp /usr/share/common-licenses/BSD .
    head BSD | sed 's/Redistribution/Mangling/'
    sed '4q;d' BSD

[sed tutorial](https://www.digitalocean.com/community/tutorials/the-basics-of-using-the-sed-stream-editor-to-manipulate-text-in-linux)

---

**4\.** for loops

    for x in /share/biocore-archive/Leveau_J_UCD/RNASeq_Arabidopsis_2016/00-RawData/*/*; do basename $x; done
    for x in /share/biocore-archive/Leveau_J_UCD/RNASeq_Arabidopsis_2016/00-RawData/*/*; do NAME=`basename $x`; echo $NAME is a file; done

---

**5\.** find

    find /share/biocore/joshi/projects/genomes -name "*.fa"
    find /share/biocore/joshi/projects/genomes -name "*.fa" -exec ls -l {} \;

---

**6\.** xargs

    find /share/biocore/joshi/projects/genomes -name "*.fa" | xargs ls -lh

---

**7\.** Installing software and git

    git clone https://github.com/najoshi/sickle.git
    cd sickle
    cat Makefile
    make

Now you need to edit your PATH.

---

**8\.** .bashrc/.bash_profile and aliases setup & the PATH variable

---

**9\.** Intro to perl

---

**10\.** Intro to mysql

    wget https://ucdavis-bioinformatics-training.github.io/2017-June-RNA-Seq-Workshop/monday/db.sqlite3
    sqlite3 db.sqlite3

---

**11\.** Intro to python


