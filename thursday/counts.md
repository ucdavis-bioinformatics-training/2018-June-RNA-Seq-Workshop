From Alignments to Raw Counts
==============================

In this section, we will collate all of the count data into one file for analysis in R.

**1\.** First, go back to your 03-alignment directory. Let's use a wildcard to list all of the counts files from all of the STAR alignment directories:

    cd ~/rnasq_example/03-alignment
    ls -1 *_star_alignment/*ReadsPerGene.out.tab

Take a look at the beginning of one of these files:

    head C61_S67_star_alignment/C61_S67_ReadsPerGene.out.tab

The first four lines are totals and the columns are ID, total reads, reads mapped to forward strand, and reads mapped to the reverse strand. In this experiment, it looks like the reads are from the reverse strand, due to the much higher mapping numbers in that column. So what we what is just that column of numbers (minus the first four lines), for every one of these files.

---

**2\.** So let's take one file and figure out how to do that, then we will expand it to all the files. First let's just get the rows we want, i.e. everything but the first four:

    tail -n +5 C61_S67_star_alignment/C61_S67_ReadsPerGene.out.tab | head

When you giv the 'tail' command a number preceded by a '+' sign, it gives you the entire file starting at the line indicated by the number. In this case, we want to skip the first 4 lines, so we start at line 5. We're piping the command to 'head' just to check that it looks correct. You shouldn't see the first four total lines.

Now, we want only the fourth column (the counts), and in order to get that we pipe the output of the tail command to the 'cut' command, and then redirect the output to a new file:

    tail -n +5 C61_S67_star_alignment/C61_S67_ReadsPerGene.out.tab | cut -f4 > C61_S67_star_alignment/C61_S67_ReadsPerGene.out.tab.count

Now, C61_S67_ReadsPerGene.out.tab.count contains a single column of data... counts for each of the genes for that sample.

---

**3\.** Now, we want to do these steps for ALL of the read count files... and to do that we will be using a 'for loop' directly on the command line. First, just run a simple 'for loop' that will print out the names of all the files we want to use:

    for x in *_star_alignment/*ReadsPerGene.out.tab; do echo $x; done

This command takes all the files that we listed in step 1 and loops through them, one by one, and for every iteration, assigns the filename to the '$x' variable. Also, for every iteration, it runs whatever commands are between the 'do' and 'done'.... and every iteration the value of '$x' changes. The semi-colons separate the parts of the loop. The 'echo' command just prints the value of $x to the screen... in this case just the filename. However, instead, we will use our previously created command, but with $x instead of the filename, and adding a few things:

    cd ../    # make sure you're in the dir above 03-alignment
    mkdir 04-Counts
    for x in 03-alignment/*/*ReadsPerGene.out.tab; do \
        s=`basename $x | cut -f1 -d_`
        echo $s
        cat $x | tail -n +5 | cut -f4 > 04-Counts/$s.count
    done

After this command, there should be a counts file for every sample in 04-Counts.

---

**4\.** Next, we need to get the columns for the final table. Because all of these files are sorted in the exact same order (by gene ID), we can just use the columns from any of the files:

    tail -n +5 03-alignment/C61_S67_star_alignment/C61_S67_ReadsPerGene.out.tab | cut -f1 > geneids.txt
    head geneids.txt

Finally, we want to combine all of these columns together using the 'paste' command, and put it in a temporary file:

    paste geneids.txt 04-Counts/*.count > tmp.out

---

**5\.** The final step is to create a header for our final counts file and combine it with the temp file. The header is just all of the sample names separated by tabs. But also, let's only take the first part of the sample name, because we actually don't need the second part. And again, since we pasted the columns in sorted order (wildcards automatically sort in order), the columns just need to be in that same order... which is the order in our samples.txt file.

    cat samples.txt | cut -d_ -f1 | paste -s > header.txt

NOTE: Workaround for working with a symbolically linked 03-alignment directory:

    for x in 03-alignment/*/*ReadsPerGene.out.tab; do \
        s=`basename $x | cut -f1 -d_`
        echo $s
    done | paste -s > header.txt


We take the samples.txt file, cut out the first column where the delimiter is the underscore character, then pipe that to the 'paste' command with the '-s' option, which takes a column of values and transposes them into a row. And finally, let's put everything together:

    cat header.txt tmp.out > all_counts.txt
    
And now you have a raw counts file that has a count for every gene, per sample. You will use this file for the next step, which is analysis in R.
