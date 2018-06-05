Preprocessing Data
===================

In this exercise, we will learn how to preprocess our data for alignment. We will be doing adapter trimming and quality trimming. Make sure you are logged into a compute node, not the head node (cabernet).

**1\.** First, create a directory for the example in your home directory:

    cd
    mkdir rnaseq_example

---

**2\.** Next, go into that directory and link to the directory for your raw data. This data comes from an Arabidopsis RNA-Seq project that we did:

    cd rnaseq_example
    ln -s /share/biocore-archive/Leveau_J_UCD/RNASeq_Arabidopsis_2016/00-RawData

---

**3\.** Now, take a look inside that directory.

    cd 00-RawData
    ls

--- 

**4\.** You will see a list of directories and some other files. Take a look at all the files in all of the directories:

    ls *

---

**5\.** Pick a directory and go into it. Look at one of the files using the 'zless' command:

    cd I894_S90_L006
    zless I894_S90_L006_R1_001.fastq.gz

Make sure you can identify which lines correspond to a single read and which lines are the header, sequence, and quality values. Press 'q' to exit this screen. Then, let's figure out the number of reads in this file. A simple way to do that is to count the number of lines and divide by 4 (because the record of each read uses 4 lines). In order to do this, use "zcat" to output the uncompressed file and pipe that to "wc" to count the number of lines:

    zcat I894_S90_L006_R1_001.fastq.gz | wc -l

Divide this number by 4 and you have the number of reads in this file. One more thing to try is to figure out the length of the reads without counting each nucleotide. First get the first 4 lines of the file (i.e. the first record):

    zcat I894_S90_L006_R1_001.fastq.gz | head -4

Then, copy and paste the sequence line into the following command (replace [sequence] with the line):

    echo -n [sequence] | wc -c

This will give you the length of the read. See if you can figure out how this command works.

---

**6\.** Now go back to your 'rnaseq_example' directory and create another directory called '01-Trimming':

    cd ~/rnaseq_example
    mkdir 01-Trimming

---

**7\.** Go into that directory (make sure your prompt shows that you are in the 01-Trimming directory) and link to all of the files you will be using:

    cd 01-Trimming
    ln -s ../00-RawData/*/*.gz .
    ls -l

Now you should see a long listing of all the links you just created.

---

**8\.** Before we do any trimming, let's run a quality control check on one of the files. To do this, we will use a piece of software called 'FastQC'. Load the module and check out the usage & options:

    module load fastqc
    fastqc -h

FastQC creates html output that has graphics with quality control analysis. You'll need to create an output directory first and then run fastqc:

    mkdir fastqc_out
    fastqc -t 6 -o fastqc_out I894_S90_L006_R1_001.fastq.gz

When that is done, you will need to download the fastqc_out directory to your laptop and then use a browser look at the html file in it. We will go over the contents of the output in class.

---

**9\.** Now, we will use software called 'scythe' (developed at the UC Davis Bioinformatics Core) to do adapter trimming. First we will run it on just one pair of files. First, load the module, and then type 'scythe' with no arguments to see the options.

    module load scythe
    scythe

Looking at the Usage you can see that scythe needs an adapter file and the sequence file. The adapter file will depend upon which kit you used... typically you can find the adapters from the sequencing provider. In this case, Illumina TruSeq adapters were used, so we have put the adapters (forward & reverse) in a file for you already ([adapters file](adapters.fasta)). You will have to use the "wget" command to copy the file to your class directory:

    wget https://ucdavis-bioinformatics-training.github.io/2017-June-RNA-Seq-Workshop/tuesday/adapters.fasta

Now run scythe specifying an output file, the adapters file, and the input file. Add an ampersand at the end to run it in the background so that we can run the other file through scythe concurrently:

    scythe -o I894_S90.scythe.R1.fastq -a adapters.fasta I894_S90_L006_R1_001.fastq.gz &
    scythe -o I894_S90.scythe.R2.fastq -a adapters.fasta I894_S90_L006_R2_001.fastq.gz &

This will take approximately 5 minutes to run. You can use the 'top' or 'jobs' commands to monitor the jobs. When the jobs finish, you will have two files that are adapter trimmed.

---

**10\.** Once that is done, let's take a look at the differences between the input and output files. First look at the input file:

    zless I894_S90_L006_R1_001.fastq.gz

Let's search for the adapter sequence. Type '/' (a forward slash), and then type **AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC** (the first part of the forward adapter). Press Enter. This will search for the sequence in the file and highlight each time it is found. You can now type "n" to cycle through the places where it is found. When you are done, type "q" to exit. Now look at the output file:

    less I894_S90.scythe.R1.fastq

If you scroll through the data (using the spacebar), you will see that some of the sequences have been trimmed. Now, try searching for **AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC** again. You shouldn't find it. You may need to use Control-C to get out of the search and then "q" to exit the 'less' screen.

---

**11\.** Now we will trim for low quality using a program called 'sickle' (also developed at the Core). First, load the module and type 'sickle' by itself to get the usage, and then 'sickle pe' to get the usage for paired-end reads.

    module load sickle
    sickle
    sickle pe

Our reads are paired-end reads in separate files, so we will be using the "-f", "-r", "-o", and "-p" options. Remember that you will be using the scythe output files as input to this step.

    sickle pe -f I894_S90.scythe.R1.fastq -r I894_S90.scythe.R2.fastq -o I894_S90.sickle.R1.fastq -p I894_S90.sickle.R2.fastq -s I894_S90.singles.fastq -t sanger

This will take about 5 minutes to run. If you look through the output files, you will see reads trimmed for low quality. Sickle produces three files, two paired-end quality trimmed files, and a singles file where reads are kept where only one of the pair passed the thresholds. Sickle will output information about how many records it started with and how many were kept/discarded.

---

**12\.** We have run through adapter & quality trimming for one pair of files, but now we need to do it for all the files. For that we will be using the power of our cluster. You'll need to logout of your compute node and get back to the head node (cabernet). You'll need to download two cluster scripts [qa_task_array.sh](qa_task_array.sh) and [qa_for_loop.sh](qa_for_loop.sh) :

    wget https://ucdavis-bioinformatics-training.github.io/2017-June-RNA-Seq-Workshop/tuesday/qa_task_array.sh
    wget https://ucdavis-bioinformatics-training.github.io/2017-June-RNA-Seq-Workshop/tuesday/qa_for_loop.sh
    
We will also need to generate a file called 'samples.txt' that contains all the sample IDs. We will extract this information from the filenames using 'cut'. First, get a listing of all the R1 files:

    ls -1 *R1*.fastq.gz

We will pipe this output to 'cut' to get the fields we want. Give cut the options "-d" with an underscore (usually above the minus sign on a keyboard) as the parameter, and the "-f" option with "1,2" as the parameter in order to get the first and second fields:

    ls -1 *R1*.fastq.gz | cut -d_ -f1,2

This gives us the all the sample IDs. Now we just need to redirect that output to a file:

    ls -1 *R1*.fastq.gz | cut -d_ -f1,2 > samples.txt

Use 'cat' to view the contents of the file to make sure it looks right:

    cat samples.txt

---

**13\.** There are many different ways to run jobs on the cluster and on the command-line... we are going to talk about two of the ways. Let's take a look at the two scripts we downloaded. The first is a script that uses Slurm task arrays to run all of the sickle and scythe steps per sample. The second is a script that uses a 'for loop' to loop through all of the samples and run the steps serially. This second script can be used when you are running all of your jobs on one machine. Look at the first script:

    cat qa_task_array.sh

You will see that it has a few extra sbatch options. The main option to understand is the "--array" option. This option creates a "task array" to run jobs. What that means is that Slurm will run this job however many times specified (in this case, 24) and for every time it runs this script will assign an environment variable called "$SLURM_ARRAY_TASK_ID". This variable will get assigned the number 1 for the first time the script runs, the number 2 the second time the script runs, etc... all the way to 24 (we are using 24 because there are 24 samples). This number is then used as an index into the samples.txt file that you created earlier. The command used to get the sample name is 'sed'. 'sed' is a program that does text editing, but in this case we are using it to get the Nth line of the samples.txt file. So, for example, this command:

    sed "5q;d" samples.txt

will return the 5th line of the samples.txt file. We put the command in backticks (usually below the tilde on a keyboard) which tells the script to run the command and put the output into the 'sample' variable. And instead of "5", we use the $SLURM_ARRAY_TASK_ID variable that will change for every run of the script. So, in effect, what happens is that the script gets run 24 times and each time the $SLURM_ARRAY_TASK_ID variable is assigned a new number, which is then used to get the sample ID from the samples.txt file.

---

**14\.** Take a look at the other script:

    cat qa_for_loop.sh

This script has similar commands, but instead of using a task array, it is using a for loop. So this will loop through all the IDs in samples.txt and assign a new ID on every iteration of the loop. You should use this script if you will be running jobs NOT on a cluster, but on a single machine.

---

**15\.** However, since we ARE running on a cluster we will use the first script to run all our jobs. Now, this step will take many hours to run, so you should probably only run it at the end of the day. First, make sure the script is executable:

    chmod a+x qa_task_array.sh
    
Now, since we have already set up the task array to run, all we need to do is run the script:

    sbatch qa_task_array.sh

Now you can use 'squeue' to make sure your jobs are queued properly. Now, all you have to do is wait.
