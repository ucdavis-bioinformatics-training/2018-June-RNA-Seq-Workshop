Advanced Command-Line
=======================

The sed command
----------------

Let's take a closer look at the 'sed' command. sed (short for stream editor) is a command that allows you to manipulate character data in various ways. One useful thing it can do is substitution. First, make a directory called "advanced" in your home directory and go into it.

    cd
    mkdir advanced
    cd advanced

Let's copy over a simple file to work on:

    cp /usr/share/common-licenses/BSD .

Take a look at the file:

    cat BSD

Now, let's change every occurence of the word "Redistribution" into "Mangling":

    cat BSD | sed 's/Redistribution/Mangling/gi'

Let's break down the argument to sed (within the single quotes)... The "s" means "substitute", the word between the 1st and 2nd forward slashes (i.e. /) is the word the substitute for, the word between the 2nd and 3rd slashes is the word to substitute with, and finally the "gi" at the end are flags for global substitution (i.e. substituting along an entire line instead of just the first occurence on a line), and for case insenstivity (i.e. it will ignore the case of the letters when doing the substitution).

Note that this **doesn't** change the file itself, it is simply piping the output of the cat command to sed and outputting to the screen. If you wanted to change the file itself, you could use the "-i" option to sed:

    cat BSD
    sed -i 's/Redistribution/Mangling/gi' BSD

Now if you look at the file, the lines have changed.

    cat BSD

Another useful use of sed is for capturing certain lines from a file. You can select certain lines from a file:

    sed '4q;d' BSD

This will just select the 4th line from the file.

**CHALLENGE:**
See if you can find a way to use sed to remove all the spaces from the BSD file.

More pipes
-----------

Now, let's delve into pipes a little more. Pipes are a very powerful way to look at and manipulate complex data using a series of simple programs. First take a look at the contents of the "/home" directory:

    ls /home

These are all the home directories on the system. Now let's say we wanted to find out how many directory names begin with each letter. First we cut out the first letter of the directories:

    ls /home | cut -c1

In order to do the counting, we first need to sort the data and then send it to the "uniq" command to keep only the unique occurences of a letter. The "-c" option counts up the number of occurences:

    ls /home | cut -c1 | sort | uniq -c

You'll notice there might be a letter missing... why would that be?

Now let's look at some fastq files. Link a few files into the advanced directory:

    ln -s /share/biocore-archive/Leveau_J_UCD/RNASeq_Arabidopsis_2016/00-RawData/C61_S67_L006/C61_S67_L006_R*_001.fastq.gz .
    
Since the files are gzipped files we need to use "zcat" to look at them. zcat is just like cat except for gzipped files:

    zcat C61_S67_L006_R1_001.fastq.gz | head

Notice that each header line has the barcode for that read at the end of the line. Let's count the number of each barcode. In order to do that we need to just capture the header lines from this file. We can use "sed" to do that:

    zcat C61_S67_L006_R1_001.fastq.gz | sed -n '1~4p' | head

By default sed prints every line. In this case we are giving the "-n" option to sed which will **not** print every line. Instead, we are giving it the argument "1~4p", which means to print the first line, then skip 4 lines and print again, and then continue to do that.

Now that we have a way to get just the headers, we need to isolate the part of the header that is the barcode. There are multiple ways to do this... we will use the cut command:

    zcat C61_S67_L006_R1_001.fastq.gz | sed -n '1~4p' | cut -d: -f10 | head

So we are using the "-d" option to cut with ":" as the argument to that option, meaning that we will be using the delimiter ":" to split the input. Then we use the "-f" option with argument "10", meaning that we want the 10th field after the split. In this case, that is the barcode.

Finally, as before, we need to sort the data and then use "uniq -c" to count. Then put it all together and run it on the entire dataset (This will take about a minute to run):

    zcat C61_S67_L006_R1_001.fastq.gz | sed -n '1~4p' | cut -d: -f10 | sort | uniq -c

Now you have a list of how many reads were categorized into each barcode. Here is a [sed tutorial](https://www.digitalocean.com/community/tutorials/the-basics-of-using-the-sed-stream-editor-to-manipulate-text-in-linux) for more exercises.

One final thing to know is that if a program does not take input from STDIN (which is needed to use it in a pipe), but instead wants a filename, you can use a single dash by itself in place of the filename and the shell will interpret that to be input from STDIN. So it would look something like this:

    cat FILENAME | COMMAND -f - -otheroptions | ....

**CHALLENGE:**
Find the distribution of the first 5 bases of all the reads in C61_S67_L006_R1_001.fastq.gz. I.e., count the number of times the first 5 bases of every read occurs across all reads.

Process substitution
---------------------

Next, we will cover process substitution. Process substitution is a way of using the output of some software as the input file to another software without having to create intermediate files. Let's use the sickle program we compiled earlier (and should already be in our PATH). We want to do adapter trimming on one of our fastq.gz files, but we need to give sickle an uncompressed file as input. In order to do that, we use the "gunzip" command with the "-c" option. This unzips the file and sends the output to STDOUT, instead of unzipping the file in place which is the default (This will take a few minutes to run):

    sickle se -f <(gunzip -c C61_S67_L006_R1_001.fastq.gz) -t sanger -o trimmed.fa

So we are putting the gunzip command inside parentheses with a less-than symbol like so: <(COMMAND). When we do this, the output of the COMMAND gets manipulated by the shell so that sickle thinks it is a file. Sickle then uses this "file" as the input file. Take a look at the output file:

    less trimmed.fa

Loops
------

Loops are useful for quickly telling the shell to perform one operation after another, in series. For example:

    for i in {1..21}; do echo $i >> a; done  # put multiple lines of code on one line, each line terminated by ';'
    cat a
    # <1 through 21 on separate lines>

The general form is:

    for name in {list}; do
        commands
    done

The list can be a sequence of numbers or letters, or a group of files specified with wildcard characters:

    for i in {3,2,1,liftoff}; do echo $i; done  # needs more excitement!
    for i in {3,2,1,"liftoff!"}; do echo $i; done  # exclamation point will confuse the shell unless quoted
    # Now imagine you have 20 sequence files, in a 'fastqs' directory:
    bwa index reference.fa
    for sample in fastqs/*.fastq; do
        bwa mem reference.fa $sample 1> $sample.sam 2> $sample.err
    done
    # this would produce, for example, ./fastqs/sample1.fastq.sam and ./fastqs/sample1.fastq.err, etc.

Sometimes a "while" loop is more convenient than a "for" loop ... if you don't readily know how many iterations of the loop you want:

    while {condition}; do
        commands
    done

Or, imagining a file that contains the filenames (one per line) of samples' sequence data:

    cat file-of-filenames.txt | while read sample; do
        bwa mem reference.fa $sample 1> $sample.sam 2> $sample.err
    done

Now, let's use a for loop on some fastq files. You can specify all the parts of the for loop on one line, separated by semi-colons. First let's just echo all the names of the fastq files in a directory:

    for x in /share/biocore-archive/Leveau_J_UCD/RNASeq_Arabidopsis_2016/00-RawData/*/*; do echo $x; done

Now, let's use the "basename" command to get just the filename for each file:

    for x in /share/biocore-archive/Leveau_J_UCD/RNASeq_Arabidopsis_2016/00-RawData/*/*; do basename $x; done

We can also assign the output of a command to a new variable by using the backtick character (\`). We put the command we want inside backticks and then assign it to the variable NAME. Then we can use $NAME in the next command:

    for x in /share/biocore-archive/Leveau_J_UCD/RNASeq_Arabidopsis_2016/00-RawData/*/*; do NAME=`basename $x`; echo $NAME is a file; done

We can also use the backticks to generate the list for the for loop. Let's say we wanted to iterate over just the sample names in /share/biocore-archive/Leveau_J_UCD/RNASeq_Arabidopsis_2016/00-RawData. First, let's generate a list of the sample names and put them in a file:

    ls -d /share/biocore-archive/Leveau_J_UCD/RNASeq_Arabidopsis_2016/00-RawData/*_L006 | cut -f7 -d/ > samples.txt

This command gets the directories ending in "\_L006" and then cuts out the 7th field using the "/" as the delimiter. Take a look at the file:

    cat samples.txt

Now, we will use this file to generate the list in the for loop by using the backticks:

    for x in `cat samples.txt`; do echo Do something with $x; done

**HARD CHALLENGE:**
Use a for loop with pipes to recreate the result from above where we wanted to find how many directory names in /home began with each letter. You will need to create a for loop to get the letters and then pipe the result of the for loop to commands to do the counting.

Find
-----

Find is a very powerful command that is used to recursively find files/directories in a file system. First take a look at the find man page:

    man find

Notice there are LOTS of options. The simplest version of the find command will simply list every file and directory within a path, as far down as it can go:

    find /share/biocore/joshi/projects/genomes

If we want to refine the command to only show you files that end in ".fa" (i.e. fasta files), we use the "-name" option:

    find /share/biocore/joshi/projects/genomes -name "*.fa"

One of the most powerful uses of find is to execute commands on every file it finds. To do this, you use the "-exec" option. When you use that option, everything after the "-exec" is assumed to be a command, and you use the "{}" characters to substitute for the file names that it finds. So in the command below, "wc -l" will get executed sequentially for every file it finds. Finally, the exec option needs to end with a semi-colon, however, since the semi-colon is a special character that the shell will try to interpret, you need to "escape" the semi-colon with a backslash, to indicate to the shell that the semi-colon is NOT to be interpreted and just sent as is to the find command:

    find /share/biocore/joshi/projects/genomes -name "*.fa" -exec wc -l {} \;

You will probably want to Ctrl-C out of this because it will take a long time to go through them all.

Xargs
------

xargs is another command that can be very useful for running a program on a long list of files. For example, the "find" commands we ran above could be run using xargs like this:

    find /share/biocore/joshi/projects/genomes -name "*.fa" | xargs wc -l

This is taking the output of the find command and then creating a list of all the filenames which it adds to the command given to xargs. So, in this case, after "xargs" comes "wc -l"... so "wc -l" will get run on the entire list of filenames from find.


.bashrc/.bash_profile, aliases & the PATH variable
-----------------------------------------------------

On a Linux system, there is usually a user-modifiable file of commands that gets run every time you log in. This is used to set up your environment the way that you want it. On our systems, the file is ".bash_profile" and it resides in your home directory. Sometimes the file is called ".bashrc" as well. Take a look at a .bash_profile:

    cat /home/joshi/.bash_profile

This one has a lot of stuff in it to set up the environment. Now take a look at your own .bash_profile. One of the things you set up was adding to the PATH variable. One thing that is very useful to add to the PATH variable is the "." directory. This allows you to execute things that are in your current directory, specified by ".". So, let's use nano to edit our .bash_profile and add "." to PATH:

    nano ~/.bash_profile

Add ":." to the PATH variable. Now, next time you log in, "." will be in your PATH.

Another thing that is very useful are aliases. An alias is a user-defined command that is a shortcut for another command. For example, let's say you typed the command "ls -ltrh" a lot and it would be easier to have it be a simpler command. Use an alias:

    alias lt='ls -ltrh'

Now, you've created an alias that lists the contents of a directory with extra information (-l), in reverse time order of last modified (-r and -t), and with human readable file sizes (-h). Try it out:

    lt

Typing alias by itself give you a list of all the aliases:

    alias

You can now put the alias command for lt in your .bash_profile and you will have it automatically when you log in.

More grep
----------

Grep is a very powerful tool that has many applications. Grep can be used to find lines in a file that match a pattern, but also you can get lines above and below the matching line as well. The "-A" option to grep is used to specify number of lines after a match and the "-B" option is for number of lines before a match. So, for example, if you wanted to find a particular sequence in a fastq file and also the 3 other lines that form that entire fastq record, you would do this:

    zcat C61_S67_L006_R1_001.fastq.gz | grep -B1 -A2 CACAATGTTTCTGCTGCCTGAACC

This looks for the sequence "CACAATGTTTCTGCTGCCTGAACC" in the fastq file and then also prints the line before and two lines after each match. 

Another thing grep can do is regular expressions. Regular expressions are a way of specifying a search pattern. Two very useful characters in regular expressions are "^" and "$". The "^" symbol specifies the beginning of a line and the "$" specifies the end of a line. So, for example, if you wanted to find just the lines that began with "TTCCAACACA" you would do this:

    zcat C61_S67_L006_R1_001.fastq.gz | grep ^TTCCAACACA

Without the "^", grep will find any line that has "TTCCAACACA" *anywhere* in the line, not just the beginning. Conversely, if you wanted to find the lines that ended in "TAAACTTA":

    zcat C61_S67_L006_R1_001.fastq.gz | grep TAAACTTA$

There are also extended regular expression that grep can use to do more complex matches, using the "-E" option:

    zcat C61_S67_L006_R1_001.fastq.gz | grep -E '^TTCCAACACA|TAAACTTA$'

This command will find any line that begins with "TTCCAACACA" **OR** ends with "TAAACTTA". The "\|" character means OR.

**CHALLENGE:**
Find a way to use grep to match any line that has between 7 and 16 'A's in a row at the end of the line. You will probably need to look at the man page for grep.

The Prompt
------------

The Prompt is the part of the command line that is the information on the screen before where you type commands. Typically, it has your username, the name of the machine you are on, and your current directory. This, like everything else in Linux, is highly customizable. Your current prompt probably has your username, your hostname, and your current directory. However, there are many other things you could add to it if you so desired. The environment variable used for your prompt is "PS1". So to change your prompt you just need to reassign PS1. There are some special escape characters that you use to specify your username, hostname, etc... these can be found in the "PROMPTING" section of the bash man page:

    man bash

Then type "/" and "PROMPTING" to find the section. You'll see that for username it is "\u", for hostname it is "\h", and for the full current working directory it is "\w". So you would change your prompt by doing this:

    export PS1="\u@\h:\w\\$ "

The convention is to put the "@" symbol and the ":" symbol to delineate the different parts, but you can use anything you want. Also, it is convention to put a "$" at the end to indicate the end of the prompt.

You can also do all kinds of fancy things in your prompt, like color and highlighting. Here is an example of such a prompt:

    export PS1='\[\033[7;29m\]\u@\h\[\033[0m\]:\[\e[1m\]\w\[\e[m\]$ '

When you have a prompt you like, you can put it in your .bash_profile/.bashrc so that it is automatically set when you log in.

Nohup
------

The nohup (short for "no hangup") command is useful for running a job from a terminal and then wanting to exit the terminal. When you run a job, even if you put it in the background (i.e. by using "&"), the job is tied to the terminal you are running on. When you log out of that terminal, any job tied to that terminal will be killed. This is not desirable, so you can use the nohup command which will disconnect the proccess from the terminal. Simply put "nohup" in front of the command, and you will probably want to add the "&" at the end to put it in the background so you can get your prompt back. It would look something like this:

    nohup YOUR COMMAND &

Awk
----

Awk is a simple programming language that can be used to do more complex filtering of data. Awk has many capabilities, and we are only going to touch on one of them here. One really useful thing is to filter lines of a file based on the value in a column. Let's get a file with some data:

    wget https://ucdavis-bioinformatics-training.github.io/2018-March-Bioinformatics-Prerequisites/tuesday/DMR.GBM2.vs.NB1.bed

Take a look at the beginning of the file:

    head DMR.GBM2.vs.NB1.bed

Let's say we wanted to get only the lines where the pvalue (column 10) was below a certain value. In awk, the default delimiter is the tab character, which most bioinformatics files use as their delimiter. Using the tab as a delimiter, it assigns each value in a column to the variables $1, $2, $3, etc... for as many columns as there are. So if we wanted to get lines where the pvalue column was under 0.00000005 pvalue:

    cat DMR.GBM2.vs.NB1.bed | awk '$10 < 0.00000005'

And lines where the pvalue >= 0.00000005:

    cat DMR.GBM2.vs.NB1.bed | awk '$10 >= 0.00000005'

You can also use it to extract lines with a particular word in a column:

    cat DMR.GBM2.vs.NB1.bed | awk '$1 == "chr3"'

A double equals (==) is used for equality comparisons. This will pull out lines where the chromosome column is "chr3".

Take a look at the [awk manual](https://www.gnu.org/software/gawk/manual/gawk.html) to learn more about the capabilities of awk.

**HARD CHALLENGE**:
Go through the list of genomes (as in the Find section) and this time only search down a maximum of 6 directories and also follow symbolic links in the search. Then extract only those files that are part of either the zebrafish or C. elegans genomes. For each of those files, get the number of characters in the file and then only print files whose character count is greater than 10000. You will have to probably use find, grep, xargs, wc, and awk. You will need to look at the manual pages for each of those commands. You should be able to do this just using pipes and the commands (i.e. no intermediate files).

