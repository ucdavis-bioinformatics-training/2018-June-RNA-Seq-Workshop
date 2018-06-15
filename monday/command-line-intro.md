Intro to the Command Line
==========================

Getting There
--------------------------------

Secure SHell ... SSH. Replace [class#@] with your login name in the following:

    ssh [class#@]ganesh.genomecenter.ucdavis.edu
    # for example:
    ssh class42@ganesh.genomecenter.ucdavis.edu

To start off: there will be many commands that will fill your screen with text. There are multiple ways to clear the clutter, and have an empty screen:

    <enter> <enter> <enter> 
    <ctrl-l>
    <ctrl-k? for Macs>
    clear

And once you're really done working on the command line:

    exit  # kills the current shell!
    # Note that any text *following* a '#' symbol is ignored.
    # This will 'exit' (close) your shell on ganesh
    # ... you'd need another 'exit' to close your local terminal (if on a Mac / Linux system)

But don\'t exit yet ... or if you did, just ssh back in ... We\'ve got work to do!

Command Line Basics
--------------------

First some basics - how to look at your surroundings.

    pwd  # present working directory ... where am I?
    ls   # list files here ... you should see nothing since your 'class##' homes are empty
    ls /tmp/  # list files somewhere else

Let\'s run our first command ... because one of the first things that\'s good to know is *how to escape once you\'ve started something you don\'t want*.

    sleep 1000  # wait for 1000 seconds!
    <ctrl-c>  # shows as '^C'; exits command entry without executing
    # In some cases, a different key sequence is required
    python  # enter Python interactive session
    <ctrl-d>  # escape from some repl\'s ... 'Read, Execute, & Print Loop'
    R  # enter R interactive session
    <ctrl-d>
    # Try this fun command:
    curl -s http://metaphorpsum.com/sentences/1000  # that's some deep stuff
    curl -s http://metaphorpsum.com/sentences/1000 | more  # 'pipe (|),' or direct, the text to the 'more' command
    <spacebar>  # next page
    <q>  # quits from more, less, 'man' pages, etc.

So, ^C, ^D, 'q', and (from above) 'exit'. Generally can't hurt to try until one of them works! Now that we've seen the 'less' *paginator* (program that display output from other programs in pages, instead of all at once), here's how to move around while in it.

    curl -s http://metaphorpsum.com/sentences/1000 | fold -w 80 -s | less  # pipe to 'less' paginator, instead of 'more'
    <spacebar>  # next page
    <arrow keys, pgup, pgdn>  # forward or back through file
    g  # beginning ...
    G  # ... or end of file
    /A<enter>  # '/' enters search mode, 'A' is pattern looked for (could be any string)
    /;<enter>  # find semicolons
    /.<enter>  # surprise! '.' is a special character that means any character
    /T.e<enter>  # should find 'The', 'They', maybe 'Toe'?
    # After searching for some text, using the '/' key, you can use:
    n  # next pattern match ...
    N  # ... or previous pattern match
    q  # to quit


You've Got Options
-------------------

One reason you'll appreciate 'less' is that it's the default paginator for the 'man' command. 'man' stands for 'manual', and it's the main way to get more detail on any of the commands we'll introduce today. Each command can act as a basic tool, or you can add 'options' or 'flags' that modify the default behavior of the tool. These flags come in the form of '-v' ... or, when it's a more descriptive word, two dashes: '\-\-verbose' ... that's a common (but not universal) one that tells a tool that you want it to give you output with more detail. Sometimes, options require specifying amounts or strings, like '-o results.txt' or '\-\-output results.txt' ... or '-n 4' or '\-\-numCPUs 4'. Let's try some, and see what the man page for the 'list files' command 'ls' is like.

    ls -R /software
    # ack! too much going to the screen!
    <ctrl-c>
    ls -R /software/scythe  # lists directories and files *recursively*
    # how do I know which options do what?
    man ls
    # navigate like in 'less' (up,down,pgup,pgdn,g,G,/pattern,n,N,q)
    # look up and try the following:
    ls -l  # if you don't say where, it lists files in your current directory
    ls -a
    ls -l -a
    ls -la  # option 'smushing' ... when no values need specifying
    ls -ltrha
    ls -ltrha --color  # single letter (smushed) vs word options
    # Quick aside: what if I want to use same options repeatedly? and be lazy?
    alias  # lists current *command aliases* or *shortcuts*
    alias ll='ls -ltrhaF --color'
    ll
    alias l='ls --color'
    l


Getting Around
----------------

The filesystem you're working on is like the branching root system of a tree. The top level, right at the root of the tree, is called the 'root' directory, specified by '/' ... which is the divider for directory addresses, or 'paths'. We move around using the 'change directory' command, 'cd':

    cd  # no effect? that's because by itself it sends you home (to ~)
    cd /  # go to root of tree's root system
    cd home  # go to where everyone's homes are
    pwd
    cd class42  # use your actual home, not class42
    pwd
    cd /
    pwd
    cd ~  # a shortcut to home, from anywhere
    pwd
    cd .  # '.' always means *this* directory
    pwd
    cd ..  # '..' always means *one directory up*
    pwd

Don't get confused between the '.' directory name and filenames that start with the '.' character. The latter are just valid filenames or directory names, but are usually hidden (use ls's '-a' option to see hidden files).

Absolute and Relative Paths
----------------------------
  
The sequence above was probably confusing, if you're not used to navigating filesystems this way. You can think of paths like addresses. You can tell your friend how to go to a particular store *from where they are currently* (a 'relative' path), or *from the main Interstate Highway that everyone uses* (in this case, the root of the filesystem, '/' ... this is an 'absolute' path). Both are valid. But absolute paths can't be confused, because they tell you where to start off, and all the steps along the way. Relative paths, on the other hand, could be totally wrong for your friend *if you assume they're somewhere they're not*. With this in mind, let's try a few more:

    cd ~  # let's start at home
    cd ../../home/class42/  # *relative* (start here, take two steps up, then down through home and class42)
    pwd
    cd /home/class42/  # *absolute* (start at root, take steps)
    pwd

Linux also tolerates 'empty' steps and loops, even if they look ugly. So:

    cd  # starting at your home again
    cd /home//class42/  # used two slashes; it's treated as an empty step
    cd ../class42/../class42/../class42/  # back and forth a few times.
    cd /software/bwa/../bowtie/../../home/class42/  # are you lost or something?

There's no real point to such weird paths, but it helps illustrate how paths work. In the last example, an absolute path, we start at the root of the filesystem (the initial '/'), move down into the 'software' directory, then down from there into the 'bwa' directory, then back up (into the software directory), down into 'bowtie', then up twice (which gets you back to the root directory), then down through the familiar 'home' and your own 'class##' directories. 

Note that the shell can be picky about how you describe a directory, or what's *in* a directory. The following two commands illustrate this difference, but just realize that you may have to experiment sometimes to get the behavior you want:

    ls -la /software  # tells you about the software directory itself
    ls -la /software/  # tells you what's *in* the software directory

Now, because it can be a real pain to type out these long paths, we need to discuss ... 

Tab Completion - A Real Tendon-Saver
-------------------------------------

Using tab-completion will literally save your life. Hours of it. A single <tab> auto-completes file or directory names when there's only one name that could be completed correctly. If multiple files could satisfy the tab-completion, then nothing will happen after the first <tab>. In this case, press <tab> a second time to list all the possible completing names. Note that if you've already made a mistake that means that no files will ever be completed correctly from its current state, then <tab>'s will do nothing.

    touch one seven september  # create three empty files using 'touch' command
    cat o<tab><no enter>  # will complete to 'one'
    <enter>
    cat s<tab><no enter>  # completes up to 'se' since that's in common between seven and september
    <tab><no enter>  # this second tab should cause listing of seven and september
    v<tab><no enter>  # now it's unique to-, and should complete to seven
    <enter>  # runs 'cat seven' command
    # we often literally autocomplete commands letter by letter
    # comes in handy if we don't exactly remember what name we want
    ls /hom<tab>j<tab><tab>f<tab>  # completes to my home directory, /home/jfass/

I can't overstate how useful tab completion is. You should get used to using it constantly. Watch experienced users type and they maniacally hit tab once or twice in between almost every character. You don't have to go that far, of course, but get used to constantly getting feedback from hitting tab and you will save yourself a huge amount of typing and trying to remember weird directory and filenames.

CHALLENGE
-----------

After returning to your home directory (just enter 'cd' by itself), verify that the two following commands are equivalent (replacing 'class42' with your actual username):

    cd ../../home/class42/; pwd  # pwd gives you your Present Working Directory
    cd ../../../../../../../home/class42/; pwd

Why are these very different-looking commands equivalent??

Create and Destroy
-------------------

OK, so let's get down to actually making some changes to the filesystem.

    cd  # home again
    mkdir temp  # make a directory called 'temp'
    cd temp/
    echo 'Hello, world!' > first.txt  # push text into a file using the '>' character
    file first.txt  # tells us what kind of file it is
    cat first.txt  # 'cat' means 'concatenate'

If a file isn't text, you probably don't want to look at it. Binary ('data') files that get pushed to the screen by the 'cat' command are chewed up in character sized bites, so your terminal displays whichever (unintended) character that bite corresponds to in the ASCII table ... which can be weird non-printing characters like bells and interrupts that can really muck up your shell! So, stick to 'cat'ing text files.
    
    # why 'concatenate'? try this:
    cat first.txt first.txt first.txt > second.txt
    cat second.txt
    # OK, let's destroy what we just created:
    cd ../
    rmdir temp  # 'rmdir' meands 'remove directory', but this shouldn't work!
    rm temp/first.txt
    rm temp/second.txt  # clear directory first
    rmdir temp  # should succeed now

So, 'mkdir' and 'rmdir' to create and destroy (empty) directories. 'rm' to remove files. To create a file can be as simple as using 'echo' and the '>' (redirection) character to put text into a file. Even simpler is the 'touch' command.

    touch newFile
    ls -ltra  # look at the time listed for the file you just created
    cat newFile  # it's empty!
    sleep 60  # go grab some coffee
    touch newFile
    ls -ltra  # same time?

So 'touch' creates empty files, or updates the 'last modified' time. Note that the options on the 'ls' command you used here give you a Long listing, of All files, in Reverse Time order (l, a, r, t).

Piping and Redirection
-----------------------

Pipes ('\|') allow commands to hand output to other commands, and redirection characters ('>' and '>>') allow you to put output into files.

    mkdir CLB
    cd CLB/
    echo 'first' > test.txt
    cat test.txt
    echo 'second' > test.txt
    cat test.txt
    echo 'third' >> test.txt
    cat test.txt

The '>' character redirects output of a command that would normally go to the screen instead into a specified file. '>' replaces, '>>' appends.

    cut -c 1-3 test.txt  # cuts character one to three, from every line, from file 'test.txt'
    cat test.txt | cut -c 1-3  # same thing, piping output of one command into input of another
    cat test.txt | cut -c 1-3 | sort -r
    cat test.txt | cut -c 1-3 | sort -r | grep s
    # pipes cat to cut to sort (-r means reverse order sort, grep searches for pattern matches)

This is a great way to build up a set of operations while inspecting the output of each step in turn. We'll do more of this in a bit.

REALLY CHALLENGING CHALLENGE
-------------------------------

We (the Bioinformatics Core) install all of the biology-related tools in the '/software/' directory, in our spare time. Are we more productive at any time of the year? So, the idea is to look at the dates the software was installed and find when the most installations occurred.

* HINT #1: Use the 'man' pages for tools we've seen already.
* HINT #2: Sometimes it's useful to turn several of the same character into only one ... for example:

    tr "LOOOOL" "LOL"

But maybe there's a way to make this behavior more general? See the man page for the translate command 'tr'.
* HINT #3: Got a sorted list? Squeeze all repeated entries into one using the 'uniq' command. But maybe 'uniq' can do more than just squeeze repeats.


History Repeats Itself
-----------------------

Linux remembers everything you've done (at least in the current shell session), which allows you to pull steps from your history, potentially modify them, and redo them. This can obviously save a lot of time and typing.

    <up arrow>  # last command
    <up>  # next-to-last command
    <down>  # last command, again
    <down>  # current command, empty or otherwise
    # try however many ups and downs you want, to step around in your history
    # for a more global view:
    history  # usually too much for one screen, so ...
    history | head
    history | tail
    history | tail -n 30
    history | less
    cat test.txt | cut -c 1-3 | sort -r | grep s > reallyImportantResult.txt
    history | tail
    !560  # re-executes 560th command (yours will have different numbers; choose your last one)

You can also search your history from the command line:

    <ctrl-r>fir  # should find most recent command containing 'fir' string: echo 'first' > test.txt
    <enter>  # to run command
    <ctrl-c>  # get out of recursive search
    <ctr-r>  # repeat <ctrl-r> to find successively older string matches

CHALLENGE
-------------

What's the first command you executed today? How many times have you used the 'man' command today? Whatever that number is, it should be more! Just kidding. Sort of.

Editing Yourself
-----------------

Here are some more ways to make editing previous commands, or novel commands that you're building up, easier:

    <up><up>  # go to some previous command
    <ctrl-a>  # go to the beginning of the line
    <ctrl-e>  # go to the end of the line
    # now use left and right to move to a single word (surrounded by whitespace)
    <ctrl-k>  # delete from here to end of line
    <ctrl-w>  # delete from here to beginning of preceeding word
    blah blah blah<ctrl-w><ctrl-w>  # leaves you with only one 'blah'


Compression
-------------

With BIG DATA(TM), you'll often see compressed files, or whole compressed folders.

    gzip test.txt
    file test.txt.gz
    gunzip -c test.txt.gz | more  # '-c' leaves the original file alone, but dumps expanded output to screen
    ls -ltrha  # see? nothing's changed
    gunzip test.txt.gz  # now the file should change
    bzip2 test.txt; bunzip2 test.txt.bz2  # note the ';' is a substitute for <enter>


Archives
--------

Tape archives, or .tar files, are one way to compress entire folders and all contained folders into one file. When they're further compressed they're called 'tarballs'. Let's grab one.

    wget ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/PhiX/Illumina/RTA/PhiX_Illumina_RTA.tar.gz
    # ... or ...
    curl ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/PhiX/Illumina/RTA/PhiX_Illumina_RTA.tar.gz > PhiX.tgz
    # .tar.gz and .tgz are *commonly used* extensions for compressed tar files, when gzip compression is used
    tar -xzvf PhiX_Illumina_RTA.tar.gz  # or whatever you called it, if you used curl
    # -x = extract, -z = use gzip/gunzip, -v = verbose (show each file in archive), -f = use a file, not a tape drive(!)

Note that, unlike Windows, linux does not depend on file extensions to determine file behavior (most of the time). So you could name a tarball 'fish.puppy' and the extract command above should work just fine. In fact, go ahead and try this. The only thing that should be different is that tab-completion doesn't work within the 'tar' command if it doesn't see the 'correct' file extension.    

Forced Removal
---------------

This gets a heading all its own. Because when you're on the command line, there's no 'Recycle Bin'. Since we've expanded a whole directory tree, we need to be able to quickly remove a directory without clearing each subdirectory and using 'rmdir'.

    rm -rf PhiX  # be sure ... there's no going back!
    # -r = recursively remove sub-directories, -f means *force* (auto-'yes')
    # We actually want to use those directories, so un-archive them again!

Obviously, be careful with 'rm -rf'. 


BASH Wildcard Characters and Find
----------------------------------

When we want to specify or operate on sets of files all at once.

    ls ?hiX/Illumina  # list files in Illumina sub-directory of any directory ending in 'hiX'
    ls PhiX/Illumina/RTA/Sequence/*/*.fa  # list all .fa files a few directories down
    # So, '?' fills in for zero or one character, '*' fills in for zero or more characters
    find . -name "*.fa"
    find . -name "*.f?"  # how is this different from the previous command?

SOMEWHAT CHALLENGING CHALLENGE
---------------------------------

Many types of software, including GNU/Linux itself, have directories named 'bin' at various levels, which are meant to hold the 'binary', compiled, executable tools themselves (as opposed to notes about the tools, data files, source code used to create the binaries, etc.). See if you can count how many 'bin' directories are in the /software directories.

* HINT #1: The 'find' man page is wild and wooly; try searching for the *second* occurrence of the text '-type c'.
* HINT #2: Acquaint yourself with the useful 'wc' tool.
* HINT #3: The /software directory is a little special. Try searching *underneath* it, whatever that's supposed to mean. 

Quick Note About the Quote(s)
-------------------------------

The quote characters " and ' are different. In general, single quotes preserve the *literal* meaning of all characters between them. On the other hand, double quotes allow the shell to see what's between them and make substitutions when appropriate. For example:

    VRBL=someText
    echo '$VRBL'
    echo "$VRBL"

However, some commands try to be 'smarter' about this behavior, so it's a little hard to predict what will happen in all cases. It's safest to experiment first when planning a command that depends on quoting ... list filenames first, instead of changing them, etc. Finally, the 'backtic' characters \` (same key - unSHIFTED - as the tielde ~) causes the shell to interpret what's between them as a command, and return the result.

    echo `$VRBL`  # tries to execute a command with the name *someText*
    newVRBL=`echo $VRBL`
    echo $vewVRBL

Delightfully confusing, eh? Just be prepared to experiment.

Manipulation of a FASTA File
-----------------------------

We just found the phiX-174 genome, so let's copy it to our current directory so we can play with it:

    cp ./PhiX/Illumina/RTA/Sequence/WholeGenomeFasta/genome.fa phix.fa
    # Note how we copied the 'genome.fa' file to a different name: 'phix.fa'
    wc -l phix.fa

We can use the 'grep' command to search for matches to patterns (more flexibly than by using less's '/' key). 'grep' comes from '**g**lobally search for a **r**egular **e**xpression and **p**rint'. 

    grep -c '>' phix.fa  # only one FASTA sequence entry, since only one header line ('>gi|somethingsomething...')
    cat phix.fa  # this may not be useful for anything larger than a virus!
    # let's look at start codon and 2 following:
    grep "ATG......" phix.fa  # '.' characters are the single-character wildcards for grep
    # use the '-o' option to **o*nly print the pattern matches, one per line
    grep -o "ATG......" phix.fa
    # use the 'cut' command with '-c' to select characters 4-6, the second codon
    grep -o "ATG......" phix.fa | cut -c4-6
    # 'sort' the second codon sequences (default order is same as ASCII table; see 'man ascii')
    grep -o "ATG......" phix.fa | cut -c4-6 | sort
    # combine successive identical sequences, but count them ('-c' option)
    grep -o "ATG......" phix.fa | cut -c4-6 | sort | uniq -c
    # and finally sort using only the 1st 'field' as a key ('-k1,1'), in reverse numeric order ('-rn')
    grep -o "ATG......" phix.fa | cut -c4-6 | sort | uniq -c | sort -rn -k1,1
    # ... which gives us the most common codons first

This may not be a particularly useful thing to do with a genomic FASTA file, but it illustrates the process by which one can build up a string of operations, using pipes, in order to ask quantitative questions about sequence content. More generally than that, this process allows one to ask questions about files and file contents and the operating system, and verify at each step that the process so far is working as expected. The command line is, in this sense, really a modular workflow management system.

CHALLENGE
----------

The commands above only find start codons on the forward strand. How would you find the most common *first codons* (after the ATG) on the reverse strand? BONUS: Can you add these in with the codons from the example above, and count them all at once?

Symbolic Links
---------------

Since copying or even moving large files (like sequence data) around your filesystem may be impractical, we can use links to reference 'distant' files without duplicating the data in the files. Symbolic links are disposable pointers that refer to other files, but behave like the referenced files in commands.

    ln -s PhiX/Illumina/RTA/Sequence/WholeGenomeFasta/genome.fa .
    ls -ltrhaF  # notice the symbolic link pointing at its target
    grep -c ">" genome.fa


Bioinformatics, At Last! ... ?
-------------------------------

OK, let's try to do some sequence alignment (similar to a BLAST search).

    # some version of BWA should already be installed, so:
    bwa mem genome.fa phix.fa > aln.sam
    # We get the following:
    # [E::bwa_idx_load] fail to locate the index files

What went wrong? We redirected output to the 'aln.sam' file using the '>' character. But there's nothing in the output file:

    cat aln.sam
    # <nothing>


STDOUT & STDERR
-----------------

Programs can write to two separate output streams, 'standard out' (STDOUT), and 'standard error' (STDERR). The former is generally for direct output of a program, while the latter is supposed to be used for reporting problems. I've seen some bioinformatics tools use STDERR to report summary statistics about the output, but this is probably bad practice. Default behavior in a lot of cases is to dump both STDOUT and STDERR to the screen, unless you specify otherwise. In order to nail down what goes where, and record it for posterity:

    bwa mem genome.fa phix.fa 1> aln.sam 2> aln.err
    # the 1st output, STDOUT, goes to 'aln.sam'
    # the 2nd output, STDERR, goes to 'aln.err'
    cat aln.sam
    # <nothing>
    cat aln.err
    # [E::bwa_idx_load] fail to locate the index files

That didn't really help us with our problem, but at least you can direct output where you want it now. Saving STDOUT is pretty routine (you want your results, yes?), but remember that explicitly saving STDERR is important in a computing cluster, since you may not directly see the 'screen' when you're running jobs. And even if you do, it could be a mix of STDERR from multiple jobs (if it's dumped to screen), so explicitly naming each STDERR file lets you diagnose what went wrong with which set of data, after everything's run or failed to run.

For the curious, our problem was that we didn't index the 'target' or 'reference' of our sequence alignment. Try this:

    bwa index genome.fa
    bwa mem genome.fa phix.fa 1> aln.sam 2> aln.err
    cat aln.err
    # see? summary of the alignment process / results
    cat aln.sam

The resulting SAM file shows a perfect match between the query sequence ('phix.fa') and the reference ('genome.fa') ... which is good, because they're the same file (see the symbolic links section above)! Beyond that, don't worry about the SAM format; we'll get into that later.


Paste Command
--------------

The paste command is useful in creating tables.

    echo 'WT1' >> b
    echo 'WT2' >> b
    echo 'control1' >> b
    echo 'control2' >> b
    # now we can number our four samples to conveniently refer to them in order
    for i in {1..4}; do echo $i >> a; done
    paste a b > c
    cat c


Running in the Background
--------------------------

Sometimes it's useful to continue working on the command line, even after you've executed a command that's going to take a while to finish. Normally this command would occupy the shell, and prevent you from typing in commands and receiving results. But we can 'put jobs in the background' so that they don't occupy your shell directly:

    sleep 1000000
    <control-z>  # shows as '^Z'
    # [1]+  Stopped                 sleep 1000000
    bg
    # [1]+ sleep 1000000 &

'^Z' first suspends the sleep command (or whatever command is running when you hit <control-z>). Then, 'bg' resumes running that command *in the background*, so that it doesn't occupy the terminal. The output of the 'bg' command tells you that you have one command running in the background. You could start more, suspend them, then resume them in the background, and query what background jobs are running or are suspended, not running:

    jobs
    # [1]+  Running                 sleep 1000000 &
    
We can also start a job in the background in one step, without having to suspend then resume it, using the '&' character at the end of the command:

    sleep 5000000 &

If we want to delete these jobs for any reason, we can kill them using the numbering that 'jobs' reveals:

    jobs
    # [1]-  Running                 sleep 1000000 &
    # [2]+  Running                 sleep 5000000 &
    kill %1
    jobs
    # [1]-  Terminated              sleep 1000000
    # [2]+  Running                 sleep 5000000 &
    kill %2
    jobs
    # [2]+  Terminated              sleep 5000000

Finally, the 'nohup' command (from 'no hangup'!) makes jobs extra resistant to lost connections or terminal problems. In other words, even jobs running in the background can be terminated if one's shell dies. 'nohup' separates the running job from the shell, so it'll keep running until it dies or is actively killed by you.

    nohup sleep 1000000 &
    # [1] 34993
    # class##@c4-0:~/CLB$ nohup: ignoring input and appending output to ‘nohup.out’
    jobs
    # [1]+  Running                 nohup sleep 1000000 &
    # output is dumped into the 'nohup.out' file unless specifically redirected in your command
    kill %1


Table of Processes
-------------------

The 'top' command prints a self-updating table of running processes and system stats. Use 'q' to exit top, 'z' to toggle better color contrast, 'M' to sort by memory use, 'P' to sort by processor use, and 'c' to toggle display of the full commands. Hit '1' to toggle display of all processors, and hit 'u' followed by typing in a username in order to only show processes (jobs) owned by that user.

![top-example](top.png)

Shell Scripts, File Permissions
--------------------------------

Often it's useful to define a whole string of commands to run on some input, so that (1) you can be sure you're running the same commands on all data, and (2) you don't have to type the same commands in over and over! Let's use the 'nano' text editor program that's pretty reliably installed on most linux systems (MacOS??).

    nano test.sh
    # nano now occupies the whole screen; see commands at the bottom
    # type/paste in the following ...
    # (note that '#!' is an interpreted command to the shell, not a comment)
    # (also note that you may have to add and delete spaces and <enter>s to get the spacing right, though that's not critical)
    #!/bin/bash
    grep -o . $1 | \
        sort | \
        uniq -c | \
        sort -rn -k1,1
    <control-o><control-x>  # to write **o**ut to test.sh, and then e**x**it nano

Note that '$1' means 'the value of the 1st argument to the shell script' ... in other words, the text that follows the shell script name when we run it (see below).

Though there are ways to run the commands in test.sh right now, it's generally useful to give yourself (and others) 'execute' permissions for test.sh, really making it a shell script. Note the characters in the first (left-most) field of the file listing:

    ll test.sh
    # -rw-r--r--  1 jfass biocore   79 Aug 19 15:05 test.sh

The first '-' becomes a 'd' if the 'file' is actually a directory. The next three characters represent **r**ead, **w**rite, and e**x**ecute permissions for the file owner (you), followed by three characters for users in the owner's group, followed by three characters for all other users. Run the 'chmod' command to change permissions for the 'test.sh' file, adding execute permissions ('+x') for the user (you) and your group ('ug'):

    chmod ug+x test.sh
    ll test.sh
    # -rwxr-xr-- 1 jfass biocore 79 Aug 19 15:05 test.sh*

OK! So let's run this script, feeding it the phiX genome. When we put the genome file 1st after the name of the script, this filename becomes variable '1', which the script can access by specifying '$1'.

    ./test.sh genome.fa
    #   1686 T
    #   1292 A
    #   1253 G
    #   1155 C
    #      1 x
    #      1 p
    #      1 i
    #      1 h
    #      1 >

The script's grep command splits out every character in the file on a separate line, then sorts them so it can count the occurrences of every unique character and show the most frequent characters first ... a quick and dirty way to get at GC content.


Installing (Simple) Software
-------------------------------

Let's install a straightforward tool ... another read aligner by the author of BWA that's intended for longer, lower accuracy reads (such as from the PacBio or Oxford Nanopore sequencers): minimap2.

    cd  # returns you to your home directory, since we've been working in the 'CLB' directory
    mkdir tools
    cd tools/
    git clone https://github.com/lh3/minimap2.git
    cd minimap2/
    make
    ./minimap2
    cd

Now you can run this tool (an executable file) by fully specifying the path to your 'tools' directory from wherever you're running the tool. In addition, you could copy, move, or put symbolic links to the tool in a common place, like /usr/bin/, which should be in everybody's path.


