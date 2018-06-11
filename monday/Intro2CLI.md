Introduction to CLI
=======================

**1\.** Find out where we are

    pwd

---

**2\.** Create a few new directories

    mkdir test
    mkdir test2
    mkdir test3

Alternatively, you could have created all three directories at once, using ***mkdir test test2 test3***

---

**3\.** Change directory

    cd

If you are in your home directory when you issue this command, you will notice that you have not changed places. That is because the command without any argument takes you to your home directory. If you were anywhere but your home directory, you will notice that you have changed to your home directory.

Now you can try to provide the command with a path, either absolute, or relative. For example:

    cd test
    pwd

You will see that you are inside the directory called ***test***.

    cd

You will see that you are back to your home directory.

    cd test
    pwd
    cd ../test2
    pwd

You will see that you went into the directory ***test***, then changed to be inside directory ***test2***.

---

**4\.** Create a file

    nano

After you give this command, you will see an interface looking like ![nano.1](./nano.1.png)

At the bottom of the terminal, you see all the commands that are available. To write any text, you just need to start typing. When you are finished, simply use the ***Exit*** command, which is to press the control key and X at the same time. You will be asked whether you want to save the file. Say ***yes***, and then give the name ***test.txt***, hit return and you will have created a new file.

---

**5\.** Other ways to create files

One may use ***cp*** command to copy a file to another:

    cp test.txt test2.txt

One can also use ***cat*** command to concatenate a few files into one file:

    cat test test2.txt > test3.txt

Here we see a new thing ***\>***. This means that the output of the command on the left (before the > ) is written in the file whoes name is given on the right: to re-direct the output of the command before it into a file. Otherwise, the output, by default, is to be the standard output (on screen).
