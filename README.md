# bash_commands

Quick reference for biologists setting up or using a VM/HPC environment with `sudo` privileges.  
Covers navigation, file management, compression, search/replace, scripting, job submission, software management, and nf-core usage.

---

## Navigating the CLI

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
  -- `ls -lhtr` = long format, human readable, sorted by time (newest last)  

---

## Creating and Managing Files

- `mkdir <name>` = make directory  
- `mkdir -p path/to/dir` = make nested directories  

- `touch <file>` = create an empty file  

- `cp <source> <destination>` = copy file  
- `cp -r <dir1> <dir2>` = copy directory recursively  

- `mv <source> <destination>` = move (or rename) file  

- `nano <file>` = simple text editor inside terminal  
- `nano -c <file>` = open file, show line numbers, create if missing  
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

## Working with Conda/Miniconda

- `conda create -n <envname> python=3.10` = create environment  
- `conda activate <envname>` = activate environment  
- `conda deactivate` = deactivate environment  
- `conda install -c bioconda samtools` = install bioinformatics tool  
- `conda env list` = list environments  
- `conda remove -n <envname> --all` = remove environment  

---

## Git and Version Control

- `git clone <repo>` = clone repository  
- `git pull` = update repository  
- `git status` = check repo status  
- `git add <file>` = stage file for commit  
- `git commit -m "message"` = commit staged changes  
- `git push` = push changes to remote repository  

---

## Automating with Loops
**Basic loop for files**  

- `for file in *.fastq.gz; do echo "Processing $file"; done`

- `for r1 in *R1*.fastq.gz; do r2=${r1/_R1/_R2}; echo "Paired reads: $r1 and $r2"; done`

- `nohup bash -c "for r1 in *R1*.fastq.gz; do r2=${r1/_R1/_R2}; echo "Paired reads: $r1 and $r2"; done" &`

---

## Writing Basic Bash Scripts
**Template: create a script.sh file that finds R2 for each R1 and run command**

```
#!/bin/bash

for r1 in *R1*.fastq.gz; do
  r2=${r1/_R1/_R2}
  sample=${r1%%_R1*}   # strip everything from _R1 onward
  echo "Processing $sample with $r1 and $r2"
  # Example: run fastqc
  fastqc $r1 $r2 -o ./qc_reports/
done
```

  -- Save the file and make it executable:

- chmod +x myscript.sh
- ./myscript.sh

---

## Managing Processes and Jobs

- `&` = run process in background
- `jobs -l` = list <PID> of background jobs
- `fg` = bring background job to foreground
- `bg` = resume job in background
- `kill <PID>` = stop process
- `kill -9 <PID>` = force kill process
- `ps aux | grep <name>` = search running processes

---

##Submitting Jobs with SLURM (if installed)
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
- `nextflow pull nf-core/<pipeline> = update pipeline`
- `nextflow -version = check installation`

---

## Useful Shortcuts

- `CTRL + C` = stop process
- `CTRL + Z` = suspend process
- `!!` = repeat last command
- `history` = show command history
- `!123` = rerun command number 123 from history
- `tab` = autocomplete file or command
- `clear` = clear terminal screen
