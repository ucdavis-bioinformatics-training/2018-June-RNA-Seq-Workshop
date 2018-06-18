Logging In & Transferring Files
================================

For Macs/Linux - Logging In
----------------------------

1. Open a Terminal (usually under Applications/Utilities on a Mac)
2. Cut and paste this into the terminal:

    ssh [class#]@ganesh.genomecenter.ucdavis.edu

where '[class#]' is replaced with your username. Press Enter.

3. Type in your password. No characters will display when you are typing. Press Enter.

For Macs/Linux - Transferring files
------------------------------------

1. Use scp (secure copy, a remote file copying program):

    scp [class#]@ganesh.genomecenter.ucdavis.edu:[full path to file] .

Replace '[class#]' with your username and replace '[full path to file]' with the full path to the file you want to transfer. Note that there is a "." at the end of the command, which is where to put the file, i.e. your current directory. You will have to type in your password.

For Windows - Logging In
-------------------------
Windows needs a terminal to connect with, you have a few options:

1. [Ubuntu for Windows 10](https://www.howtogeek.com/249966/how-to-install-and-use-the-linux-bash-shell-on-windows-10/)
2. [MobaXterm](https://mobaxterm.mobatek.net/)
3. [PuTTY](http://www.putty.org/)

Then,
1. Open up terminal (putty, mobaxterm, or ubuntu).
2. In the Host Name field, type **ganesh.genomecenter.ucdavis.edu**
3. Make sure the Connection Type is SSH.
4. Press "Open". It will ask you for your username and password.


For Windows - Transferring files
---------------------------------

1. Open up WinSCP. If you haven't installed it, get WinSCP [here](https://winscp.net/eng/download.php).
2. In the Host Name field, type **ganesh.genomecenter.ucdavis.edu**
2. Type in your username and password.
3. Make sure the File Protocol is SFTP.
4. Press "Login".
