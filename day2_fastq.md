# 🧬 Day 2 — FASTQ Files and Quality Control

In this session, we’ll explore **FASTQ files**, the fundamental data format of Next-Generation Sequencing (NGS). You’ll learn how to examine the structure of raw sequencing reads, assess data quality using **FastQC**, and perform read trimming and filtering with **fastp**.

---

### 🧠 Learning Objectives

- Understand the structure and components of FASTQ files.
- Learn to view and inspect sequencing reads using UNIX tools.
- Assess read quality using **FastQC**.
- Perform trimming and filtering using **fastp**.
- Compare raw and processed data using FastQC reports.

---
## 1️⃣ Create project directory
1. Login to IBU Bioinformatics sever (login8.hpc.binf.unibe.ch) using ssh protocol on VSCode.
2. Create a new directory called dataPreprocess in your course directory 
```shell
 mkdir -p course/dataPreProcess
```
a directory called raw_data and trimmed_data (to store cleaned data)
```shell
mkdir -p course/dataPreProcess/raw_data
mkdir -p course/dataPreProcess/clean_fastq_data
mkdir -p course/dataPreProcess/raw_fastqc_output
mkdir -p course/dataPreProcess/trimmed_fastqc_output
```
 and a scripts directory
 ```shell 
 mkdir -p course/dataPreProcess/scripts. 
```

### 📂 Project Directory Overview

By the end of **Day 2**, your working directory should look like this:

```text
.
├── course
│   ├── dataPreProcess
|   |   ├── raw_data
|   |   |   └── SRR1027171_1.fastq.gz
|   |   |   └── SRR1027171_2.fastq.gz
|   |   |── raw_fastqc_output
|   |   |   └── fastqc*zip
|   |   |   └── fastqc_report.html
|   |   ├── scripts
|   |   |   └── download.sh
│   |   ├── clean_fastq_data
│   |   |   └── SRR1027171_1.clean.fq.gz
│   |   |   └── SRR1027171_2.clean.fq.gz
│   |   |── trimmed_fastqc_output
│   |   |   └── fastqc_trimmed.html
│   |   |   └── fastqc*zip

```

## 2️⃣ Downloading Raw Data

1. Change the current working directory to dataPreprocess/data

```shell
cd course/dataPreProcess/raw_data
```
2. Download the following fastq files from NCBI short read archive (http://ncbi.nlm.nih.gov/sra)

```shell
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR102/001/SRR1027171/SRR1027171_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR102/001/SRR1027171/SRR1027171_2.fastq.gz
```
---
If wget didn't work please copy it from here 

```shell
cp /data/courses/pcourseb/fastq/SRR1027171_1.fastq.gz .
cp /data/courses/pcourseb/fastq/SRR1027171_2.fastq.gz .
```
-	The above files belong to the study- [GSE52194](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52194) 
-	 Click the above link to take you to SRA (Sequence Read Archive) where you can find more information about the source of the fastq files.

Questions: 
1. What is the name and size of the files you downloaded ?
2. Use unix less & head command to see how the header line of fastq looks like. 

```shell
less SRR1027171_1.fastq.gz | head
```
Questions: 
Does it look like the example you saw in the lecture ?


## 3️⃣ Running FastQC on Raw Reads
We will check the quality of the fastq files using the software program called fastqc. 

Fastqc documentation is available here 
https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/

We use the slurm cluster managment program to run the fastqc job. 
First we will try to run the fastqc job interactively with help of the srun command and then we will run it as a batch job.

### 3.1 Request resources using the following srun command 

First use the _srun_ command to request for resources like CPU and RAM to run the program 

```shell
srun --partition=pcourseb --cpus-per-task 4 --time=2:00:00 --mem=1G --pty /bin/bash
```

- srun                         # SLURM command to allocate and execute resources
- --partition=pcourseb         # Use the 'pcourseb' partition/queue
- --cpus-per-task 4          # Request 4 CPU cores
- --time=2:00:00             # Request 1 hour of runtime
- --mem=4G                   # Request 4 GB of memory
- --pty /bin/bash            # Start an interactive bash shell

Questions: 
Do you understand all the arguments passed to _srun_ ? 

### 3.2 Running fastqc at the shell prompt interactively 

Once you are logged into one of the compute nodes. (inspect your prompt). **Which node are logged into ?**

Load the software module in the following manner 

```shell
 module load FastQC/0.11.9-Java-11
 ```

 Launch fastqc with the two fastq files in the following manner 

 ```shell
  fastqc --extract SRR1027171_1.fastq.gz SRR1027171_2.fastq.gz --threads  4 -o ~/course/dataPreProcess/raw_fastqc_output
  ```
Questions: 
What is the module load command ? 

### 3.3 Running fastqc as a batch job using _sbatch_

The above job can also be launched using bash script on the head node. In this case _SLURM_ looks for node with required resources and launches it on the cluster. 
Please exit the interactive session before launching the job. Meaning type exit on the node your are in and get back to the master node.

```shell
exit
```
Are you back on the master node ? (check your prompt)

 **Now create the following script with VScode and save as 'run_fastqc.sh' in the ~/course/dataPreProcess/scripts directory** 

``` 
#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --mem=2G
#SBATCH --output=fastqc.out
#SBATCH --error=fastqc.err
#SBATCH --job-name=fastqc
#SBATCH --cpus-per-task=4
#SBATCH --partition=pcourseb

 module load FastQC/0.11.9-Java-11

fastqc --extract ~/course/dataPreProcess/raw_data/SRR1027171_1.fastq.gz ~/course/dataPreProcess/raw_data/SRR1027171_2.fastq.gz --threads  4 -o ~/course/dataPreProcess/raw_fastqc_output

```
Sumbit the job to the cluster 

```shell
sbatch run_fastqc.sh 
```
Questions: 
1. Is your job running ? (Hint :use _squeue_)
It should take ~5 mins for the job to finish. 
2. How many output files have been produced by fastqc ? what are the types ?

### 3.3 Viewing FastQC HTML Reports in VS Code

After FastQC completes, HTML reports will be in `~/course/dataPreProcess/raw_fastqc_output/`:
- `SRR1027171_1_fastqc.html`
- `SRR1027171_2_fastqc.html`

### Using Live Preview:

1. **Install Live Preview extension** (if not installed):
   - Open Extensions panel (Ctrl+Shift+X)
   - Search "Live Preview" by Microsoft
   - Install

2. **Navigate to output folder** in VS Code Explorer on the left pane:
```
   ~/course/dataPreProcess/raw_fastqc_output/
```

3. **Open HTML file**:
   - Right-click on `SRR1027171_1_fastqc.html`
   - Select "Show Preview" or "Open with Live Preview"
   - HTML will open in a VS Code tab

4. **Alternative - Open in external browser**:
   - Right-click HTML file
   - Select "Open with Live Server" (opens in your default browser)

The preview will show the complete FastQC report with all plots and quality metrics.

Answer the following Questions: 
1. How many reads did the files have in raw_fastqc_output ? 
2. Do you find any test failed ? 
3. What amount over-representation you see in the two files ?
4. What are these over-represented sequences ?

## 4️⃣ Filtering adapters and low quality bases
Several tools exist that can remove the adapters and filter low quality base. Examples include trimmomatic,fastx cutadapt, sickle etc. Here we will use fastp to remove the adapters and low quality bases. The fastp manual is available here 
https://github.com/OpenGene/fastp

With VS Code, create and save a bash script named fastp_clean.sh to clean the fastq files as follows:

``` 
#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --mem=2G
#SBATCH --output=qual.out
#SBATCH --error=qual.err
#SBATCH --job-name=fastp
#SBATCH --cpus-per-task=4
#SBATCH --partition=pcourseb

module add fastp/0.23.4-GCC-10.3.0

fastp -w 4 -q 15 -z 5 -l 50 -i ~/course/dataPreProcess/raw_data/SRR1027171_1.fastq.gz -I ~/course/dataPreProcess/raw_data/SRR1027171_2.fastq.gz -o ~/course/dataPreProcess/clean_fastq_data/SRR1027171_1.clean.fq.gz -O ~/course/dataPreProcess/clean_fastq_data/SRR1027171_2.clean.fq.gz
```
The arguments passed to cut-adapt were based on the following:  
- Trim low-quality ends from reads before adapter removal if quality is less than 15 (-q 15)
- Discard trimmed reads that are shorter than 50 bases after trimming (-l 50)
- compression level for gzip output (-z 5)
- set number of threads to 4 (-w 4) 

The best is to run the fastp algorithm using a job script. This way you are recording all the parameters and can be easily added to your methods in manuscripts for the cause reproducibility research.

## 5️⃣ Running FastQC on Clean Reads
- Now run fastqc on the cleaned fastq files. 
- write the output to trimmed_fastqc_output directory for the results.
- Use Live Preview to view the html output.
- Record the changes you see in the cleaned and trimmed reads

---

## 6️⃣  Adding Version Control to Your FASTQ Processing Workflow

1. Change the current working directory to "course"
```shell
cd ~/course
```
2. Edit `.gitignore`  file with VScode 

Add the following lines to exclude large FASTQ files and temporary data

```
# FASTQ files (too large for Git)
*.fastq.gz
*.fq.gz

# FastQC output
*_fastqc.html
*_fastqc.zip
```

3. Add and commit your directory/files:

```shell
# Add your new directory and the scripts
git add datapreProcess
git add .gitignore
# commit 
git commit -m "data processing commit: Add FASTQ processing scripts"
```
#### Exercise Questions:

1. Check the status of your repository:
```shell
git status
```
- What files are being tracked?
- What files are ignored?

2. View the Git history:
```shell
git log
```
3. Make the following improvements to your scripts:
- Add a README.md file describing the workflow
- Add comments explaining each parameter in run_fastqc or fastp_clean.sh
- Commit these changes separately

Open a file with VScode, add the following lines and save it as README.md
```
## Scripts for processing FASTQ files from the GSE52194 study:
- run_fastqc.sh: Quality control analysis using FastQC
- fastp_clean.sh: Adapter and quality trimming using fastp

## Usage:
1. Run FastQC on raw data
2. Clean reads using fastp
3. Run FastQC on cleaned data
```
# Add and commit README
```shell
git add README.md
git commit -m "Add workflow documentation"
```

```shell
# Make script improvements and commit
git add *.sh
git commit -m "Add detailed parameter documentation to scripts"
```

#### FOR the BRAVE 
4. Practice branching
Create a new branch for testing different quality parameters. Change directory to the root of the project directory (cd ~/course)
```shell
git checkout -b test/quality-params
#check the branch you are in 
git branch
```
Modify the fastp parameters in your script:
Open the file dataPreprocess/scripts/fastp_clean.sh with VScode.
- Change -q 15 to -q 20 and -l 50 to -l 60
- save the file and take a snapshot with git.

```shell
#Commit these changes
git add dataPreprocess/scripts/fastp_clean.sh
git commit -m "Test stringent quality thresholds"
# Push the new branch to GitHub
# Verify your remote 
git remote -v
git push -u origin test/quality-params
```
Go to your repository page and verify if the branch is upated.

###Questions
1. Why do we exclude FASTQ files from Git? What would be a better way to track large scientific datasets?
2. Look at the difference between your original and modified fastp parameters:

```shell
git checkout main
git diff test/quality-params dataPreprocess/scripts/fastp_clean.sh
```

## 🎉 Wrapping Up

You’ve learned how to explore, assess, and improve sequencing read quality using command-line tools. These skills form the foundation for downstream analyses such as read alignment and variant calling.

