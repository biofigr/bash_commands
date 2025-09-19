# Bash Commands

Quick reference for biologists setting up or using a VM/HPC environment with `sudo` privileges.  
Covers navigation, file management, compression, search/replace, scripting, job submission, software management, and nf-core usage.

---

## Basic Terminology

- **Terminal** = text-based interface to interact with the operating system (you type commands, see output)  
- **Shell** = the program that interprets your commands (e.g. `bash`, `zsh`)  
- **Command** = an instruction you type into the shell (e.g. `ls`)
- **CLI** = Command Line Interface
- **GUI** = Graphical User Interface  
- **Process** = a running program started by a command  
- **PID** = Process ID number (unique identifier for each running process)  
- **Job** = a process started from your current shell session (can be foreground or background)  
- **Signal** = a message sent to a process by the system or user (e.g. `CTRL+C` sends an *interrupt* signal, logout sends a *hangup* signal)  
- **stdout (standard output)** = normal output stream of a command  
- **stderr (standard error)** = error messages from a command  
- **stdin (standard input)** = input stream (keyboard or data fed into a command)  

---

## VM vs HPC

- **VM (Virtual Machine)**  
  - A single-user computer running in the cloud or on a local host  
  - You usually have **root / sudo privileges** (can install software, configure system)  
  - Good for training, testing pipelines, and learning Linux basics  
  - You control the environment (operating system, installed tools, data storage)  

- **HPC (High-Performance Computing cluster)**  
  - A shared system with many nodes (servers) used by multiple users  
  - You usually **do not** have root access — instead, use `modules` or `conda` to load tools  
  - Jobs are submitted to a **scheduler** (e.g. `SLURM`, `PBS`) rather than run directly in the shell  
  - Best for large datasets, long-running jobs, and high memory/CPU workloads  

- **Key differences for beginners**  
  - On a VM, you can type commands directly and run pipelines immediately  
  - On an HPC, you request resources and wait for the scheduler to run your job  
  - Both use the same Linux commands once you’re inside a shell  

---

## Navigating via CLI

- `/` = root directory  
- `~` = home directory  
- `.` = current directory  
- `..` = parent directory  

- `cd <path>` = change directory  
  -- `cd ~` = go to home  
  -- `cd ..` = go up one level  
  -- `cd -` = go back to previous directory  

- `pwd` = print working directory  

- `ls` = list contents  
  -- `ls -a` = include hidden files  
  -- `ls -lh` = long format, human readable sizes  
  -- `ls -ltr` = long format, sorted by time (oldest first)  

---

## Creating and Managing Files

- `mkdir <name>` = make directory  
- `mkdir -p path/to/dir` = make nested directories  

- `touch <file>` = create an empty file  

- `cp <source> <destination>` = copy file  
- `cp -r <dir1> <dir2>` = copy directory recursively  

- `mv <source> <destination>` = move (or rename) file  

- `nano <file>` = simple text editor inside terminal  
- `nano -c <file>` = open/create file, show line numbers  
- `vi <file>` = advanced text editor  

- `> <file>` = redirect output of a command into a new file  
- `>> <file>` = append output to an existing file  

---

## Deleting Files

**Deleted files are gone. No trash.**

- `rm <file>` = remove file  
- `rm -r <dir>` = remove directory recursively  
- `rm -rf <dir>` = force remove without asking  

---

## Viewing and Checking Files

- `cat <file>` = print whole file  
- `less <file>` = view file interactively, quit with `q`  
- `head <file>` = first 10 lines  
- `head -n 50 <file>` = first 50 lines  
- `tail <file>` = last 10 lines  
- `tail -n 50 <file>` = last 50 lines  
- `tail -f <file>` = follow live output (e.g. logs)  

- `wc -l <file>` = count lines in file  
- `wc -c <file>` = count bytes  

- `grep "pattern" <file>` = search for exact pattern  
- `grep -i "pattern" <file>` = case-insensitive search  
- `grep -c "pattern" <file>` = count matches  
- `grep -r "pattern" <dir>` = search recursively in directory  

- `zcat <file.gz>` = view a compressed file without decompressing  
- `gunzip -c <file.gz> | head` = decompress to screen and view top lines  

---

## Text Processing

These commands are often combined with pipes (`|`) to filter and summarize data.

- **`cut`** = extract columns from tab-delimited files  
  - `cut -f 1 file.tsv` = show the first column  
  - `cut -f 1,3 file.tsv` = show columns 1 and 3  

- **`sort`** = sort lines  
  - `sort file.txt` = alphabetical sort  
  - `sort -n file.txt` = numeric sort  
  - `sort -nr file.txt` = numeric sort, reverse order  

- **`uniq`** = remove or count duplicates (input must be sorted first)  
  - `uniq file.txt` = remove duplicate lines  
  - `uniq -c file.txt` = count occurrences of each unique line  

- **`grep`** (already introduced) = search for matching lines  
  - `grep -v "pattern" file.txt` = invert match (show lines *not* containing pattern)  
  - `grep -n "pattern" file.txt` = show line numbers  

- **Combining tools**  
  - `cut -f 1 genes.tsv | sort | uniq | wc -l`  
    = count unique gene IDs in the first column  
  - `grep -v "^@" aln.sam | cut -f 3 | sort | uniq -c | sort -nr`  
    = count how many reads map to each reference sequence in a SAM file  

---

## Search and Replace

- `sed 's/old/new/' <file>` = replace first occurrence of "old" with "new" per line  
- `sed 's/old/new/g' <file>` = replace all occurrences in each line  
- `sed -i 's/old/new/g' <file>` = replace in file directly (in-place edit)  

---

## Compression and Archiving

- `gzip <file.fa>` = compress (replaces file with `.gz`)  
- `gzip -c <file.fa> > file.fa.gz` = compress but keep original  
- `gzip -d <file.fa.gz>` = decompress  

- `tar -czvf output.tar.gz <dir>` = create compressed archive  
- `tar -xzvf file.tar.gz` = extract archive  

- `zip -r output.zip <dir>` = zip directory  
- `unzip <file.zip>` = unzip  

---

## Data Transfer

- `scp <file> user@server:/path/` = copy file to server  
- `scp -r <dir> user@server:/path/` = copy directory to server  
- `scp user@server:/path/file ./` = copy file from server  

- `rsync -avh <source> <destination>` = sync directories  

- `wget <url>` = download file from web  
- `curl -O <url>` = download file with same name  

---

## Standard Output, Error, and Background Jobs

*When you run a command in Bash, it produces two kinds of output:*

- **stdout (standard output)** = normal output (results, logs, text)
- **stderr (standard error)** = error messages

*By default, both are shown on screen.*

---

### Redirecting output

- `command > file.txt` = send stdout to a file (overwrite)
- `command >> file.txt` = append stdout to a file
- `command 2> errors.txt` = send stderr to a file
- `command > out.txt 2> err.txt` = split stdout and stderr into separate files
- `command &> all.txt` = send both stdout and stderr into the same file

---

## Background Jobs, Hangups, and Session Management

- **`&`** = run a command in the background, freeing the terminal  
  - Example: `fastqc *.fastq.gz &`  
  - Still tied to the terminal: if you log out, it stops  

- **Hangup (SIGHUP)** = signal sent to processes when the terminal closes (logout, SSH drop)  
  - Any background job without protection will die on logout  

- **`nohup`** = run a command immune to hangup  
  - Example: `nohup nextflow run nf-core/rnaseq -profile conda &`  
  - Output is saved in `nohup.out` unless redirected  

- **`disown`** = detach a background job from the shell so it won’t receive SIGHUP  

- **`screen`** = start a detachable session  
  - `screen` = start a session  
  - `Ctrl-A D` = detach  
  - `screen -r` = reattach  

- **`tmux`** = modern alternative to `screen`  
  - `tmux` = start a session  
  - `Ctrl-B D` = detach  
  - `tmux attach` = reattach  

## Environment and Modules

- `export PATH=$PATH:/path/to/bin` = add directory to PATH  
- `echo $PATH` = show current PATH  
- `which <command>` = check install path  
- `command -v <command>` = check if command exists  

- On HPCs with modules:  
  -- `module avail` = list available modules  
  -- `module load <tool>` = load tool  
  -- `module list` = show loaded modules  
  -- `module unload <tool>` = unload tool  

---

## Connecting to Remote Systems

- **`ssh` (Secure Shell)** = connect to another computer (VM or HPC)  
  - `ssh user@server` = connect to a server  
  - `ssh -p 2222 user@server` = connect on a non-default port  
  - `ssh -i key.pem user@server` = connect using an SSH key file  

- **First login examples**  
  - `ssh stephen@192.168.1.100` = login to a VM by IP address  
  - `ssh user@login.hpc.edu` = login to an HPC login node  

- **`exit`** = close the SSH session  

---

## Data Transfer

- **`scp` (secure copy)**  
  - `scp file.txt user@server:/path/` = copy file to server  
  - `scp -r dir/ user@server:/path/` = copy directory to server  
  - `scp user@server:/path/file.txt ./` = copy file from server  

- **`rsync`** (faster and can resume interrupted transfers)  
  - `rsync -avh file.txt user@server:/path/` = sync file to server  
  - `rsync -avh user@server:/path/ ./` = sync file from server  
  - `rsync -avh --progress largefile user@server:/path/` = show transfer progress  

- **Web downloads**  
  - `wget <url>` = download from web  
  - `curl -O <url>` = download keeping same filename  

---

## File Integrity and Checksums

Checksums let you confirm that files were not corrupted during download or transfer.

- **`md5sum`** = generate or check an MD5 hash  
  - `md5sum file.fastq.gz` = print hash for a file  
  - `md5sum file1.fastq.gz file2.fastq.gz` = check multiple files  
  - Compare the hash with the one provided by the source (should match exactly)  

- **Save checksums**  
  - `md5sum *.fastq.gz > checksums.md5` = create a list of hashes  
  - `md5sum -c checksums.md5` = verify all files against saved list  

---

## Pipes and Redirection

- **Redirection** = send output somewhere else instead of the screen  
  - `>` = redirect stdout (overwrite)  
    - Example: `ls > files.txt` = save file list to `files.txt`  
  - `>>` = redirect stdout (append)  
    - Example: `echo "new line" >> notes.txt`  
  - `2>` = redirect stderr (errors)  
    - Example: `command 2> errors.txt`  
  - `&>` = redirect both stdout and stderr  
    - Example: `command &> log.txt`  

- **Pipes (`|`)** = send output of one command directly into another  
  - `grep "pattern" <file> | wc -l` = count how many lines match pattern  
  - `grep "^@" file.fastq | wc -l` = count number of reads in a FASTQ file (each read starts `^` with `@`)  
  - `zcat reads.fastq.gz | head -n 8` = look at the first FASTQ record in a compressed file  
  - `samtools view -h aln.bam | grep -v "^@" | wc -l` = count mapped reads in BAM (ignoring headers)  
  - `cut -f 1 genes.tsv | sort | uniq | wc -l` = count unique gene IDs in first column  

- **Chaining multiple pipes** = you can connect many commands in sequence  
  - Example: `cat file.txt | grep "gene" | sort | uniq -c | sort -nr`  
    - find lines with `"gene"` → sort them → count unique entries → sort counts in descending order  

---

## Working with Conda/Miniconda

- **Anaconda** = full distribution of Python + Conda + hundreds of preinstalled packages (large install, not needed for HPC/VM work)  
- **Miniconda** = minimal installer for Conda (preferred: lightweight, you only install what you need)  
- **Conda** = environment and package manager (lets you install bioinformatics tools without root, and keep environments isolated)  
- **Mamba** = faster drop-in replacement for Conda (uses the same commands, but resolves dependencies much quicker)  
- **Why use Conda?** = avoids dependency conflicts, works without sudo, widely supported in bioinformatics (e.g. nf-core pipelines)  

- `conda create -n <envname> python=3.10` = create environment  
- `conda activate <envname>` = activate environment  
- `conda deactivate` = deactivate environment  
- `conda install -c bioconda samtools` = install bioinformatics tool  
- `conda env list` = list environments  
- `conda remove -n <envname> --all` = remove environment  

### Important Warning on Dependencies

- **Dependencies** = other software libraries a program needs in order to run  
- In bioinformatics, dependency conflicts are one of the **biggest headaches** (different tools require different versions of the same library)  
- **Only install programs with caution** — extra packages can break existing environments  
- This is why we use a **clean VM** or isolated Conda environments: so broken dependencies do not affect the entire system  

---

## Git and Version Control

- **What is Git?** = a system for tracking changes in files (especially code or scripts)  
- **Why use it?** = lets you save versions, collaborate, and roll back mistakes  

### Setup (first time)
- `git config --global user.name "Your Name"` = set your name  
- `git config --global user.email "you@example.com"` = set your email  
- `git config --list` = check current config  

### Basic Workflow
- `git clone <repo>` = download a repository from GitHub (or another remote)  
- `git pull` = update your local copy with remote changes  
- `git status` = show current changes  
- `git add <file>` = stage a file for commit  
- `git commit -m "message"` = save a snapshot of staged changes  
- `git push` = upload commits to remote repository  

### Branching (recommended)
- `git branch` = list branches  
- `git checkout -b newbranch` = create and switch to a new branch  
- `git checkout main` = switch back to main branch  
- `git merge newbranch` = merge a branch into current branch  

### Ignoring Files
- Create a file called `.gitignore` and list files or patterns to exclude from version control  
  - Example contents:  
    ```
    *.fastq.gz
    *.bam
    results/
    ```  

### Typical Cycle
1. Edit files  
2. `git add <file>`  
3. `git commit -m "short message"`  
4. `git push`  

---

## Variables and the `$` Symbol

- **Variable** = a name that stores a value (text, number, filename)  
  - Example: `sample=reads_R1.fastq.gz`  

- **`$`** = used to access the value stored in a variable  
  - `echo $sample` = prints `reads_R1.fastq.gz`  

- **In loops** = `$file`, `$r1`, `$r2` represent the current item in the loop  
  - `for file in *.fastq.gz; do echo $file; done`  
    = prints each FASTQ filename  

- **Command substitution** = run a command and use its output  
  - `DATE=$(date)` = run `date` and save result in variable `DATE`  
  - `echo $DATE` = prints the stored date/time  

- **Quoting variables** = always use `"$var"` when filenames may contain spaces  
  - `echo "$file"` (safe) vs. `echo $file` (can break if filename has spaces)  

---

## Automating with Loops

- **What is a loop?** = a way to repeat the same command on many files without typing it each time  

### General Syntax
- `for var in list; do command $var; done`  

### Examples

- Loop over FASTQ files:  
  - `for file in *.fastq.gz; do echo "Processing $file"; done`  

- Loop over paired-end reads:  
  - `for r1 in *R1*.fastq.gz; do r2=${r1/_R1/_R2}; echo "Pair: $r1 and $r2"; done`  

- Run FastQC on all FASTQs:  
  - `for file in *.fastq.gz; do fastqc $file -o qc_reports/; done`  

- Process paired reads with a program:  
  - `for r1 in *R1*.fastq.gz; do r2=${r1/_R1/_R2}; bbmap.sh in1=$r1 in2=$r2 out=${r1%%_R1*}.sam; done`  

### Background and Long Jobs
- Add `&` to run loop in the background  
- Example:  
  - `nohup bash -c 'for file in *.fastq.gz; do fastqc $file -o qc_reports/; done' &`  

### Notes
- Always test with `echo` first before running destructive commands  
- Use `{}` or quotes around variables if filenames have spaces:  
  - `echo "$file"` instead of `echo $file`  

---

## Writing Basic Bash Scripts

- **What is a script?**  
  - A plain text file containing a list of commands  
  - Lets you automate repetitive tasks instead of typing them manually  
  - Saved with `.sh` extension by convention  

- **Shebang line**  
  - First line should be: `#!/bin/bash`  
  - Tells the system to run the file using the Bash shell  

- **Comments**  
  - Lines starting with `#` are ignored by Bash  
  - Use comments to document what your script does  

### Example Script

**script.sh**  
```bash
#!/bin/bash

for r1 in *R1*.fastq.gz; do
  r2=${r1/_R1/_R2}          # find matching R2 filename
  sample=${r1%%_R1*}        # strip everything from _R1 onward
  echo "Processing $sample with $r1 and $r2"
  fastqc "$r1" "$r2" -o ./qc_reports/   # run FastQC on each pair
done
```

### Save the file and make it executable:

- `chmod +x myscript.sh`
- `./myscript.sh`

- `bash script.sh` = run script without making it executable 

---

## Managing Processes and Jobs

- `&` = run process in background
- `jobs -l` = list background jobs with job number and PID
- `fg` = bring background job to foreground
- `bg` = resume job in background
- `kill <PID>` = stop process
- `kill -9 <PID>` = force kill process
- `ps aux | grep <name>` = search running processes

---

## Submitting Jobs with SLURM (if installed)

**SLURM is a job scheduler for HPCs. You write a script (.sh) and submit it.**

  -- Example my_job.sh:
```
#!/bin/bash
#SBATCH --job-name=MyAlignment
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=my_job.out
#SBATCH --error=my_job.err

module load samtools
module load bwa

bwa mem genome.fa reads_R1.fq.gz reads_R2.fq.gz > aligned.sam
samtools view -b aligned.sam > aligned.bam
```

**Commands:**
  -- `sbatch my_job.sh` = submit job
  -- `squeue -u <username>` = check jobs
  -- `scancel <jobID>` = cancel job

---

## Running Pipelines (nf-core / Nextflow)
**nf-core is the Gold-Standard of reproducible bioinformatics pipelines**

- `nextflow run nf-core/rnaseq -profile conda`
- `nextflow run nf-core/sarek -profile docker`
- `nextflow pull nf-core/<pipeline>` = update pipeline
- `nextflow -version` = check installation

---

## Useful Shortcuts

- `CTRL + C` = stop process
- `CTRL + Z` = suspend process
- `!!` = repeat last command
- `history` = show command history
- `!123` = rerun command number 123 from history
- `tab` = autocomplete file or command
- `clear` = clear terminal screen
