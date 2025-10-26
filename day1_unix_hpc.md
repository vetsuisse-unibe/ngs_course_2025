---
layout: default
title: Day 1 — UNIX Exercises
---

# UNIX Exercises

In this practical, you will log in to the **Bioinformatics server** (`login8.hpc.binf.unibe.ch`) and perform all the exercises directly on the server.

To log in, you need to be connected to the **University of Bern** network via **VPN** or **eduroam**.

---
## eduroam Connection

If you are connected to the **eduroam** network, you do **not** need to use the VPN. Otherwise see below for VPN connection instructions.

---

## VPN Connection

Please log in to the UniBE VPN server using **FortiConnect**.  
👉 [FortiConnect VPN instructions](http://bit.ly/47pGZyL)

---

## SSH into the Bioinformatics Server

After connecting to either the VPN or eduroam network, use **Visual Studio Code’s Remote SSH extension** to access the server.

Check your login and password details in the **students-passwords** text file on the ILIAS repository:  
🔗 [ILIAS Student Passwords](https://ilias.unibe.ch/goto_ilias3_unibe_file_3222705_download.html)

### Steps to Connect

1. Start your **Visual Studio Code** application.  
2. Click on the **Remote Explorer** icon (shown below) on the left side of the window.  

   ![VSC](images/vsc-1.png)
3. A sub-window opens. Click **Connect to Host**.  

   ![VSC-2](images/vsc-2.png)
4. In the next window, click **Add New SSH Host**.  

   ![VSC-3](images/vsc-3.png)
5. Enter the login details as shown below (replace with your assigned student ID).  

   ![VSC-4](images/vsc-4.png)
6. Click **Connect** and enter your password when prompted.  

   ![VSC-5](images/vsc-5.png)
7. Open a new terminal by selecting **Terminal → New Terminal**.  

   ![VSC-6](images/vsc-6.png)

---

### Notes

💡 **Ask for assistance** at any point if something is unclear — the goal is to learn by doing!  

- Commands shown in grey blocks should be typed in the **Unix terminal**.  
- Words in *italics* or enclosed in `<angle brackets>` (e.g., `<filename>`) should be replaced by your own parameters such as file names or paths.

---

## 1️⃣ Try Some Basic UNIX Commands

Display your user name  
whoami

Show the current working directory  
pwd

It should show something like:  
/home/student36  
This is your **home directory**.

List all files and directories  
ls

List with details  
ls -l

List all files including hidden ones  
ls -a

---

## 2️⃣ Navigating the File System

Change directory  
cd <directory_name>

Go up one directory  
cd ..

Return to your home directory  
cd

Show your current location again  
pwd

---

## 3️⃣ Working with Files and Directories

Create a new directory  
mkdir test_dir

Change into it  
cd test_dir

Create an empty file  
touch myfile.txt

View the file  
cat myfile.txt

Copy the file  
cp myfile.txt copy.txt

Rename or move the file  
mv copy.txt renamed.txt

Delete the file  
rm renamed.txt

Go up one level and remove the directory  
cd ..  
rmdir test_dir

---

## 4️⃣ Viewing and Searching File Contents

Display the first 10 lines  
head <filename>

Display the last 10 lines  
tail <filename>

Search for a word  
grep "pattern" <filename>

Count lines, words, and characters  
wc <filename>

Display a file page by page  
less <filename>  
(Press `q` to quit)

---

## 5️⃣ Redirection and Pipes

Redirect output to a file  
ls -l > list.txt

Append output to a file  
echo "New line" >> list.txt

View the output file  
cat list.txt

Combine commands with pipes  
ls -l | grep ".txt"

Count how many `.txt` files exist  
ls -1 *.txt | wc -l

---

## 6️⃣ File Permissions

Check file permissions  
ls -l

Change permissions (add execute permission)  
chmod +x script.sh

Remove write permission  
chmod -w myfile.txt

Make a file readable and writable only to you  
chmod 600 private.txt

---

## 7️⃣ Compression and Archiving

Compress a file  
gzip myfile.txt

Decompress it  
gunzip myfile.txt.gz

Create a tar archive  
tar -cvf myarchive.tar directory_name/

Extract a tar archive  
tar -xvf myarchive.tar

Compress a tar archive  
tar -czvf myarchive.tar.gz directory_name/

Extract a compressed tar archive  
tar -xzvf myarchive.tar.gz

---

## 8️⃣ Disk Usage and File Size

Check your current disk usage  
df -h

Check how much space a folder takes  
du -sh <folder_name>

Sort files by size  
ls -lhS

---

## 9️⃣ Bonus: Useful Commands

history – Show command history  
man <command> – Show manual page for a command  
clear – Clear the terminal screen  
date – Show the current date and time  
uptime – Show system uptime  
top – Monitor processes in real time  
df -h – Show disk usage  
du -sh <dir> – Show directory size

---

## ✅ Summary

In this exercise, you learned how to:

- Connect to the Bioinformatics server  
- Navigate and manipulate files in UNIX  
- View and search file contents  
- Use redirection, pipes, and permissions  
- Compress, extract, and inspect file sizes  

These are the **essential UNIX skills** for working on the HPC cluster and performing NGS data analysis in later sessions.
```