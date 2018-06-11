Introduction to CLI
=======================

**1\.** Find out where we are

    pwd

---

**2\.** Create a few new directories

First, we are going to create a directory, where we are going to work through this tutorial.

    mkdir Intro2CLI
    cd Intro2CLI

Once we are inside the directory "Intro2CLI", we are going to practice more in creating new directories.

    mkdir test
    mkdir test2

Alternatively, you could have created both directories at once, using ***mkdir test test2***

---

**3\.** Get files from a remote location

    wget https://github.com/ucdavis-bioinformatics-training/2018-June-RNA-Seq-Workshop/blob/master/monday/art.part1.txt
    wget https://github.com/ucdavis-bioinformatics-training/2018-June-RNA-Seq-Workshop/blob/master/monday/art.part2.txt
    wget https://github.com/ucdavis-bioinformatics-training/2018-June-RNA-Seq-Workshop/blob/master/monday/art.part3.txt
    wget https://raw.github.com/ucdavis-bioinformatics-training/2018-June-RNA-Seq-Workshop/master/monday/all_counts.txt

These are three files that we are going to use later in the tutorial.

---

**4\.** Change directory

    cd

If you are in your home directory when you issue this command, you will notice that you have not changed places. That is because the command without any argument takes you to your home directory. If you were anywhere but your home directory, you will notice that you have changed to your home directory.

Now you can try to provide the command with a path, either absolute, or relative. For example:

    cd ./Intro2CLI/test
    pwd

You will see that you are inside the directory called ***test***, which is inside the directory ***Intro2CLI***, which in turn is in your home directory.

    cd

You will see that you are back to your home directory.

    cd ./Intro2CLI/test
    pwd
    cd ../test2
    pwd

You will see that you went into the directory ***test***, then changed to be inside directory ***test2***, both of these two directories are inside the directory ***Intro2CLI***.

---

**5\.** Create a file

First, we change our current location into the directory ***Intro2CLI***.

    cd ../

Then, type the following command:

    nano

After you give this command, you will see an interface looking like ![nano.1](./nano.1.png)

At the bottom of the terminal, you see all the commands that are available. To write any text, you just need to start typing. When you are finished, simply use the ***Exit*** command, which is to press the control key and X at the same time. You will be asked whether you want to save the file. Say ***yes***, and then give the name ***test.txt***, hit return and you will have created a new file.

This is one of the text editors that you may use to edit a text file in command. There are many others, such as vi/vim, emacs, ...

---

**6\.** Other ways to create files

One may use ***cp*** command to copy a file to another:

    cp test.txt test2.txt

One can also use ***cat*** command to concatenate a few files into one file:

    cat art.part1.txt art.part2.txt art.part3.txt > oliver.art.txt

Here we see a new thing: **>**. This means that the output of the command on the left (before the > ) is written in the file whoes name is given on the right: to re-direct the output of the command before it into a file. Otherwise, the output, by default, is to be the standard output (on screen). For example, if we use the command below, it simply show the content of the file on screen.

    cat oliver.art.txt

Since we have sine one way to redirect the output of a command to a file, there is another related method, by using **>>**. The difference between these two ways is that: **>** takes the standard ouput and write into a file. If the file already exist, it will overwrite the original content. While **>>** will take the standard output and append it into the file.

Now that we have learned how to create a file, the next step is to learn a few commands that will alow us to manipulate a file.

**7\.** Commands to look at a file

First, in many cases, we could like to take a look at just the first a few lines of a big files. To open the file using a text editor is always one way to do it. But when the size of a file gets bigger, it takes time to open a file using a text editor. The way to do it fast is to use the command ***head***.

    head -n 10 oliver.art.txt

The option **-n** followd by a number (N) tells the command to show the first N lines of a file.

As one can imagine, there is a command to look at the last number of lines of a file, which is ***tail***.

    tail -n 10 oliver.art.txt

There is one very useful command that can show a file page by page, which does not have to read the entire file before displaying, therefore it has advantage when looking at a hugh file.

    less all_counts.txt

While stay inside the command, we can search for any pattern that we are interested in, such as a specific gene (ATCG01090). We can do it by issue the command ***/ATCG01090***, followed by hitting the enter key on your keyboard. In order to get out of looking at the file, we simply type **q**. This will give us back the prompt and we can then issue a new command.

---

**8\.** Search for a pattern in a file

In many cases, we want to find all the appearances of a certain pattern. The command we can use is grep.

    grep 'ATCG01090' all_counts.txt

The result of this command is all the lines that match the pattern.

In the case where we have a few patterns that we are interested, we can use an option in grep ***-E***.

    grep -E "ATCG01090|AT1G03997" all_counts.txt

The option ***-E*** tells ***grep*** to use an Extended Regular Expression, where **|** is an logical operator of OR, meaning that we want to search any line that matches the pattern "ATCG01090", or "AT1G03997". The result of this command is the two lines that each matching one pattern. The command ***egrep*** is exactly the same as ***grep -E***.

This command is very good at finding matches for a few patterns. However, sometimes, we might have many more patterns that we want to search. Under this situation, we can use a different option ***-F*** and provide a file that have all the patterns that we are interested in.

    grep -F -f pattern.txt all_counts.txt

The result of this command is all the lines that match any of the patterns in our list in the file pattern.txt. The command ***fgrep*** is the same as ***grep -F***.

**9\.** Extract specific fields from a file

In the field of bioinformatics, we have to frequently extract specific columns from a file that has a delimitor to separate the columns. We can easily achieve the goal by using the command ***cut***.

    cut -f2,4:10 all_counts.txt

The command above extracts the column 2,4 to 10 from the file "all_counts.txt". By default, the command ***cut*** uses tab as the delimitor. If the file is formated using a different delimitor, we can add the option of **-d** to specify the specific delimitor of the file. For example, for a file that uses comma as the delimitor, one would add the option of ***-d','***.

   cut -d',' -f2,4:10 all_counts.csv



