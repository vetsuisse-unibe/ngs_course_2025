# 🧬 HPC Exercises

In this session, you’ll learn how to use the **High-Performance Computing (HPC)** cluster and the **SLURM job scheduler** to run bioinformatics analyses efficiently.

---

## 1️⃣ Connecting to the HPC Cluster

> **Learn how to log in to the cluster, access resources, and navigate the environment.**

- Login into the Bioinformatics server (login8.hpc.binf.unibe.ch)

- Use the remote login extension on Visual studio code and login into the bioinformatics server with your chosen username and password.

Open the Terminal and start the exercises.
---

## 2️⃣ Exploring the HPC Environment

> *Understand the filesystem, available modules, and environment variables.*

The HPC system has many compute nodes (machines) where jobs are executed. To get an overview of their current usage, you can run:
```shell
clusterstate.sh
```
This script displays the current state of all nodes managed by the SLURM scheduler.

The output of the command is a table with the following columns:

| Column | Description |
|---------|-------------|
| `NODELIST` | Node name (machine identifier) |
| `CPUS(A/I/O/T)` | CPUs **Allocated / Idle / Other / Total** |
| `CPU_LOAD` | Average CPU load on that node |
| `FREEMEMORY` | Approximate available memory (in MB) |
| `STATE` | Current node status — e.g. `idle`, `mixed`, `down`, `completing` |


| State | Meaning |
|--------|----------|
| `idle` | Node is available for new jobs |
| `mixed` | Node is running some jobs but still has free CPUs |
| `down` | Node is offline or under maintenance |
| `completing` | Node is finishing current jobs |
| `alloc` | Node is fully allocated to running jobs |

**Questions:**
- Identify which nodes are **idle** or **down**.  
- Can you spot which nodes have the **highest CPU load** or **lowest memory**?  
- **Discuss:** Why do some nodes show as *mixed* instead of *idle* or *alloc*?

## 3️⃣ Configuring Git

> *Learn how to configure Git for version control.*

Before starting the SLURM exercises, we'll set up a _git repository_ to track our work.
git (version 2.43.7) is already installed on the HPC system. Check the version of git installed on your system.

```shell
git -v
```
**We'll use git command-line tools to create the repository.**

On the command line, git commands follow this structure: `git verb options`.  
- The `verb`  tells git what action to perform (like `commit` or `push`), 
- while `options` provide extra details to modify the command's behavior.

 Lets configure the git environment on command line. You should have to do configure this only once on any given computer. We will use _git config_ command to set these variables.These configuration variables are put into a .gitconfig file and stored:

- on Mac in /Users/\<YourMacusername\>/.gitconfig
- on Windows in C:\Users\\<YourWindowsUsername\>\AppData\Local\Git\config/.gitconfig
- on Linux in /home/\<Yourusername\>/.gticonfig

```shell
# Set your git username and email
git config --global user.name "Your username"
git config --global user.email "your.email@example.com"
```
Please use your own name and email address in the place holders. This is crucial because Git permanently embeds this information into each commit you make. In the next lesson, we will be interacting with GitHub and so the email address used should be the same as the one we will use to set up your GitHub account.

You can see where the .gitconfig is stored using the below command

```shell
git config --list --show-origin
```
## 4️⃣ Create and initialize the repository
first create a directory called course and navigate into it. in

```shell
cd
mkdir course
cd course
```
initialize the repository. This will create a .git directory in the current directory. and keep track of all the changes you make to the files in the repository. 
```shell
git init
```
create a directory called hpc-exercises and navigate into it. This is where we will be working on the exercises.
```shell
mkdir -p hpc-exercises/scripts
cd hpc-exercises
```
Create a .gitignore file to exclude non-essential HPC files from being tracked or backed up.

1. Open a new text file in VS Code (File Menu $\rightarrow$ New Text File).
2. Add the necessary exclusion lines (provided below/elsewhere).
3. Save this file as .gitignore directly within your hpc-exercises directory
```
*.out
*.err
*.txt
slurm-*.out
```
Track and save the changes to the .gitignore file in your repository. Which means we to need add and commit the file.

```shell
git add .gitignore
git commit -m "Initial commit: Add .gitignore for HPC output files"
```

You home directory  directory structure will look like this at the end of the day:
```text
/home/<student_name>
└── course
    └── hpc-exercises
        └── scripts
```
You can use the tree command to view the directory structure.
```shell
tree .
```
## 5️⃣ Submitting Batch Jobs

> *Learn how to write and submit SLURM batch scripts for reproducible analyses.*

### Exercise 1
In this exercise, we’ll write a simple Bash script using the VS Code command, run it on our local machine, and then submit it to the SLURM scheduler using the sbatch command.

1. Open a new text file in VS Code (File Menu $\rightarrow$ New Text File).
2. Copy the lines below and paste it into the file. 

```shell
#!/bin/bash
hostname
date
sleep 30
date
```
**Notes:** 
- `#!/bin/bash` is the shebang line that tells the system this is a bash script
- `hostname` prints the name of the host (node) on which the script is running
- `date` prints the current date and time
- `sleep 30` pauses the script for 30 seconds
- `date` prints the current date and time again
3. Save the file as **test.sh** under the courses/hpc-exercises/script scripts directory 
4. Make it executable using the _chmod_ command on the terminal.

```shell
ls -l test.sh
chmod +x test.sh
ls -l test.sh
```
Now that you created a new script, add it to your git repository and commit it. Doing this will help you keep track of the changes you make to the script.
```shell
git add test.sh
git commit -m "Add initial test script"
```
Run locally on the login node (for demo purposes only. Real work should be submitted to Slurm to run on computing nodes.)
```shell
./test.sh
```

Submit to Slurm
```shell
sbatch -p pcourseb  -t 00:05:00 test.sh
```
- `sbatch` is the SLURM command to submit jobs
- `-p pcourseb` is the partition to submit the job to
- `-t 00:05:00` is the time limit for the job

Check the job status
```shell
squeue -u $USER
```
- `squeue` is the SLURM command to check the status of jobs
- `-u $USER` is a standard Unix environment variable (system wide placeholder) that automatically substitutes your current login username (student43, student44, etc).

Questions:
1. What is the job ID of the submitted script.
2. Where is the output of the job ?
3. Check the _git status_ - are any new files created that aren't tracked? (use git status)
4. What does the command _git log --oneline_ do ?

```shell
git status
git log --oneline
```

### Exercise 2: Random Number Generation

In this exercise, we will create a bash script that generates random numbers and submit it to the SLURM scheduler using sbatch. We’ll use the VScode command to write the script.

1. Open a new text file in VS Code (File Menu $\rightarrow$ New Text File).
2. Copy the lines below and paste it into the file.

```shells
#!/bin/bash
for i in {1..1000}; do echo $RANDOM >>randomNumbers.txt; done
sort -n randomNumbers.txt
```
3. Save the file as **randomNumbers.sh** under the scripts directory and make it executable using the chmod command.

```shell
ls -l randomNumbers.sh
chmod +x randomNumbers.sh
ls -l randomNumbers.sh
```
Next, stage the new script file, and then commit it to your Git history.
```shell
git add randomNumbers.sh
git commit -m "Add random number generation script"
```

Submit to Slurm using _sbatch_
```shell
sbatch -p pcourseb -N 1 -n 1 --mem 100 -t 2:00:00 -o randomNumbers.out -e randomNumbers.err randomNumbers.sh
```

Check the job status
```shell
squeue -u $USER
```
**Notes:**
```text
#SBATCH -p pcourseb            # Specifies which partition/queue to run the job on
                               # In this case, it's using the 'pcourseb' partition

#SBATCH -N 1                   # Requests 1 node for this job
                               # A node is a complete computer in the cluster

#SBATCH -n 1                   # Requests 1 CPU core/task
                               # This defines how many parallel processes to run

#SBATCH --mem 8G               # Requests 8 gigabytes of RAM for the job
                               # This is the total memory allocation

#SBATCH -t 0-2:00             # Sets the time limit for the job
                               # Format is D-HH:MM (0 days, 2 hours, 0 minutes)

#SBATCH -o randomIntegers.out          # Specifies where to write standard output (stdout)
                               # Will create a file named 'randomIntegers.out'

#SBATCH -e randomIntegers.err          # Specifies where to write standard error (stderr)
                               # Will create a file named 'randomIntegers.err'
```

Questions:
1. What is the job ID of the submitted script.
2. How many files has the job submission created ?
3. Where is the output of the job ?
4. Check git ls-files to see files tracked - which files are untracked? Why?

### Exercise 3: Resource Allocation in Script
A better approach is defining resource allocation inside the shell script. This way you will not need to remember for the next time and simply re-run the analysis if required.
Create a script with embedded SLURM parameters:

Open a text file in VSCode, add the below lines and save as **randomIntegers.slurm** under the courses/hpc-exercises/script scripts directory

```shell
#!/bin/bash
#SBATCH -p pcourseb 
#SBATCH -N 1 
#SBATCH -n 1 
#SBATCH --mem 8G 
#SBATCH -t 0-2:00 
#SBATCH -o randomIntegers.out
#SBATCH -e randomIntegers.err


for i in {1..100}; do echo $RANDOM >> randomIntegers.txt; done
sort -n randomIntegers.txt
```

Track the file in git
```shell
git add randomIntegers.slurm
git commit -m "Add random number generation script with embedded SLURM parameters"
```

Submit to Slurm using _sbatch_
```shell
sbatch randomIntegers.slurm
```

Check the job status
```shell
squeue -u $USER
```

Questions:
1. What is the job ID of the submitted script.
2. what is the job name
3. on which compute node the job is running

## 5️⃣ Monitoring and Managing Jobs

> *Discover how to view job status, cancel jobs, and check resource usage.*


### Exercise 4: Job Control
Create a long-running script to practice job cancellation.

The scancel command can be used to cancel a job after its submitted. 
1. Write the following code and save the file as **test4.sh** under the courses/hpc-exercises/script scripts directory. 
2. Submit  job to SLURM. 
3. Wait for the job to start running (status R), then cancel it prematurely using the _scancel_ command.

```shell
#!/bin/bash
#SBATCH -p courseb # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of cores
#SBATCH --mem 10G # memory pool for all cores
#SBATCH -t 0-2:00 # time (D-HH:MM)
#SBATCH -o test4.out # STDOUT
#SBATCH -e test5.err # STDERR

hostname
date
sleep 600
date
```

Track the file in git
```shell
git add test4.sh
git commit -m "Add long-running script for job cancellation practice"
```

Submit to Slurm using _sbatch_
```shell
sbatch test4.sh
```

Check the job status
```shell
squeue -u $USER
```
Wait for job to start, then:
Cancel the job using _scancel_
```shell
scancel <job_id>
```
You can also use _scancel -u $USER_ to cancel all jobs submitted by the user.

### Exercise 5: Job Monitoring
Use sacct to analyze job performance

The sacct command can tell you information about both running jobs and finished jobs. It communicates with SLURM’s database of job information and can tell you lots of useful statistics about your jobs, such as how much memory and CPU they used.

When you run sacct without any arguments, it will display a summary of all completed jobs in the system. This summary may include information such as job IDs, user names, job status, start and end times, and other job-related details.

One can use the -j/–jobs flag, where it takes the job ID as the input.

Trying running sacct without parameters or with -j flag and answer the following questions:

1. How many Jobs completed
2. How much of memory did you request and how much was used ? 
3. Can we use this to reduce the amount of memory requested next time ?
4. Review the Git log - how many commits have you made ?
```shell
git log --oneline
```
## 6️⃣ Best Practices

1. Always commit your scripts before running them
2. Use meaningful commit messages
3. Don’t track output files in Git
4. Keep your scripts organized in directories
5. Document significant changes in commit messages

We will continue with more SLURM jobs and track them with git in the rest of our exercises.


---

