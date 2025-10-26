## Unix Exercises
In this Practical we will login into the Bioinformatics server (login8.hpc.binf.unibe.ch) and do all the exercises on the server.

In order to login into the Bioinformatics cluser you need to be connected to unibe network via VPN or eduroam.
### VPN
Please login into the unibe VPN server using FortiConnect.
(https://tinyurl.com/yxj8ovpd)

### eduroam
If you are connected to the eduroam network then you don't need to connect to VPN server.

### SSH into bionformatics server
After connecting to either the VPN server or eduroam network, use Visual Studio Code's Remote SSH extension to access the server.

Check the login and password details on students-passwords text file on the Ilias respository (https://ilias.unibe.ch/goto_ilias3_unibe_file_3222705_download.html).

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

### Notes

💡 **Ask for assistance** at any point if something is unclear — the goal is to learn by doing!  

- Commands shown in grey blocks should be typed in the **Unix terminal**.  
- Words in *italics* or enclosed in `<angle brackets>` (e.g., `<filename>`) should be replaced by your own parameters such as file names or paths.

---

## 1️⃣ Try Some Basic UNIX Commands
*Display user name*
```shell
whoami
```

*Show the current working directory*
```shell
pwd
```

It should show something like /home/student36. Which is 'home' directory for user student36

*List the files in the directory*
```shell
ls
```

*Create a empty file*
```shell
touch <filename>
```
replace <filename> with your own filename (of choice)

*Print a string to the screen*
```shell
echo "Hello world"
```

*Print the current date*
```shell
date
```

*See a history of all the last commands you tried*
```shell
history
```

*Get local help page of a command*
```shell
man ls
```

*Run the following commands one after another*
```shell
touch <exampleFile>
ls
```
replace <exampleFile> with your own filename (of choice)

*Print the contents of the file*
```shell
cat <exampleFile>
```
replace <exampleFile> from the before command.
It should show no lines as it is an empty file.

*Rename a file*
```shell
mv <exampleFile> <exampleFile2>
```
replace <exampleFile1> and <exampleFile2> with your own filename (of choice)

*Delete a file*
```shell
rm <exampleFile2>
```

*Create a new folder/directory*
```shell
mkdir <exampleDirectory>
```

*Create a file under the new folder/directory*
```shell
touch exampleDirectory/exampleFile
```

*Delete the folder/directory*
```shell
rmdir exampleDirectory
```

Worked? No! 😅
That’s because rmdir only removes empty directories.

*Delete the file inside the directory first*
```shell
rm exampleDirectory/exampleFile
```

*Delete the directory now*
```shell
rmdir exampleDirectory
```

# Creating and moving around Directories
*Create a directory called 'Documents' and change current directory to Documents*
```shell
mkdir Documents
cd Documents
pwd
```

*Going up one directory*
```shell
cd ..
```

*and then type to see what has happened*
```shell
pwd
```

*Go up by two directories*
```shell
cd ../..
```

*Go to home directory*
```shell
cd
```

Always type “pwd” to locate yourself

# Command arguments

Most programs in UNIX accept arguments that modify the program’s behavior. For example
List the files in longer format
```shell
ls -l
touch exampleFile1
touch exampleFile2

ls
ls -l
```

Different example parameters used with ls
| Command | Description |
|----------|-------------|
| `ls -a` | List all files, including hidden ones |
| `ls -h` | Show file sizes in human-readable format (KB, MB, GB) |
| `ls -l` | Long format — permissions, owner, size, date |
| `ls -S` | Sort files by size |
| `ls -t` | Sort files by modification time |
| `ls -1` | Display one file per line |

*Parameters for remove command*
```shell
rm exampleFile1
rm -i exampleFile2
```

Aliases are short forms used for commands.
```shell
alias rm=”rm -i”
touch <exampleFile>
rm <exampleFile>
```
So better to alias rm as rm -i to be on the safer side ?.

## 5️⃣ Redirection and Pipes.

> *Understand input/output redirection and how to connect commands with pipes.*

All the above commands sent the output if any to the screen. Instead of outputting on the screen redirection helps you to put into a file
```shell
echo "My first line" >testFile.txt
cat testFile.txt
echo "My second line" >>testFile.txt
cat testFile.txt
ls / >> ListRootDir.txt
cat ListRootDir.txt
```

# Wildcards

wildcard is a symbol that is used to represent one or more characters. Example wildcards are as follows
* Zero or more characters
? Any single character.

In the following examples you can test some of these wildcard characters

Create a new folder

```bash
mkdir wildCardTesting
```

Change directory
```shell
cd wildCardTesting
touch test1.txt
touch test2.txt
```

Create several files with a single command
```shell
touch test3.txt test4.txt test1.csv test2.csv test3.csv test4.csv
```

Count the number of files using the pipes ( | symbol)

```bash
wc –l counts the number of lines in the input
```

```shell
ls -l | wc -l
```

Try using the following wild cards
```shell
ls *
ls test*
ls *.txt
ls *.csv
ls test[1-2].txt
ls test[!3].*
```

# Subset a file with grep and awk
### 🧠 Learning Outcomes

By the end of this exercise, you will be able to:

- Use `grep`, `awk`, `cut`, and `sort` to extract and process text data.  
- Work with **tab-delimited** files such as the NCBI RefSeq summary.  
- Filter and count genome assemblies for a specific organism.  
- Combine commands using **pipes (`|`)** for efficient data processing.  
- Apply basic UNIX text manipulation to real biological datasets.  

Create a new folder called GenomeStats to run the exercise.

```shell
mkdir GenomeStats
cd GenomeStats
```

Download the text file showing the available refseq genomes for different species at NCBI https://www.ncbi.nlm.nih.gov/genome/browse/
RefSeq Genomes on NCBI are high-quality reference genome sequences that serve as standard representations for the genetics of specific organisms
```shell
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt
```

Use less to have a quick view of the file
```shell
less assembly_summary_refseq.txt
#Press q to exit the less interface and come back to prompt.
```

It is a tab delimited text file with several columns.
The second line shows the different column headers
```shell
head -n 2 assembly_summary_refseq.txt | tail -n 1
```

For further ease of reading the column headers. (Do you understand what we are doing here ?)
```shell
head -n 2 assembly_summary_refseq.txt | tail -n 1 | tr '\t' '\n' |less
#Press q to exit the less interface and come back to prompt
head -n 2 assembly_summary_refseq.txt | tail -n 1 | tr '\t' '\n' |nl |less
#Press q to exit the less interface and come back to prompt
```

We want to see how many dog assemblies have been submitted
```shell
awk -F'\t' '$8 == "Canis lupus familiaris" && $6 == "9615"' assembly_summary_refseq.txt
```

Want to keep the header line ?
```shell
awk -F'\t' 'NR == 2 || ($1 == "Canis lupus familiaris" && $6 == "9615")' assembly_summary_refseq.txt

```

- NR == 2: NR is the current line number. NR == 2 allows the header line (second line) to be printed.
- `||` : Logical OR, so it matches either the header line or lines that meet the specified condition.

Lets do some statistics on available Genomes
How many Animal and plant genomes are available
cut command in unix can be used to select columns from tab de-limited files
Only the column Group can selected using cut
```shell
cut -f25 assembly_summary_refseq.txt | less
```

Now the pipes can be used to see the number of different refseq genomes available at NCBI
```shell
cut -f25 assembly_summary_refseq.txt | sort | uniq -c

```

We can skip the header line and before counting the genome groups in the column 25
```shell
tail -n +3 assembly_summary_refseq.txt | cut -f25 | sort | uniq -c
```

- tail -n +3 assembly_summary_refseq.txt: Outputs the file starting from the third line, effectively skipping the first two lines.
- cut -f25: Extracts the 25th column from the remaining lines.
- sort: Sorts the output to prepare for counting unique values.
- uniq -c: Counts occurrences of each unique value in the 25th column

*So how many mammalian genomes are available ?*
Now use cut and pipe symbol to find the number of mammalian genomes available at NCBI. (Hint: check column 25 for vertebrate_mammalian)
cut can be used to select more columns
```shell
cut -f 8,25,28 assembly_summary_refseq.txt |less
```

which Mammalian genome has the highest GC content
```shell
cut -f 8,25,28 assembly_summary_refseq.txt | grep vertebrate_mammalian |sort -t$'\t' -nrk3 | less
cut -f 8,25,28 assembly_summary_refseq.txt | grep vertebrate_mammalian |sort -t$'\t' -nrk3 | head -n 1
```

- sort: Sorts the filtered results.
- -t$'\t': Sets the field delimiter as a tab ($'\t' is a literal tab in Unix shell).
- -nrk3: Specifies sorting by the 3rd column in numerical order (-n), in reverse order (-r), based on column 3 (-k3).
- Result: Sorts entries in descending order based on GC% values.

* Which Mammalian genome has the least GC content ? *

# Question of the day
Invertebrates generally have lower GC content compared to vertebrates.
You can check if this is true in your dataset
```shell
awk -F'\t' '$25 ~ /vertebrate_/ { total++; if ($28+0 < 40) count++ } END { print (count/total)*100, count, total }' assembly_summary_refseq.txt
awk -F'\t' '$25 ~ /invertebrate/ { total++; if ($28+0 < 40) count++ } END { print (count/total)*100, count, total }' assembly_summary_refseq.txt
```

- total++ If column 25 contains "vertebrate_" count total vertebrate genomes
- $28+0 converts column 28 (GC%) to a number.
- count++ if GC content (column 28) is < 40% Count these low-GC genomes
- Output shows three numbers:Percentage of vertebrate genomes with GC<40%, Count of low-GC genomes, Total vertebrate genomes

# Sequence Content
Create a new directory and copy the chromosome  human chromosome 22 from here: /data/courses/courseB/UnixEx
```shell
mkdir <directoryname>
cd <directoryname>
cp /data/courses/courseB/UnixEx/chr22.fa.gz  .
```

How many nucleotides are found in the entire chromosome?
```shell
less chr22.fa.gz | grep -v ">" | wc | awk '{print $3-$1}'
```

Do you understand the above command. If not try man wc for help.

How many As, Cs, Gs, Ts and Ns are found in the entire chromosome?
```shell
less chr22.fa.gz |grep -v ">" | grep -o [actgnACTGN] | sort | uniq -c
```

Search for  EcoR1 (GAATTC) site in the chr22 file
```shell
less chr22.fa.gz | grep -v ">" | grep --color "GAATTC"
```

Now can you count the number of EcoR1 site in the sequence ?

**For the brave.** Calculate the %GC content in the entire chromosome.

```shell
less chr22.fa.gz | awk '!/^>/{gc+=gsub(/[gGcC]/,""); at+=gsub(/[aAtT]/,"");} END{ printf "%.2f%%\n", (gc*100)/(gc+at) }'
```

### 🎉 *Welcome to the exciting world of Data Analysis!*