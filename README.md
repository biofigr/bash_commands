# bash_commands

A quick reference list of common bash commands for biologists starting to use the command line on a VM or HPC environment.  
Covers navigation, file management, system checks, package/software installation, and process handling.

---

## Navigating the CLI

- `/` = root directory  
- `~` = home directory  
- `.` = current directory  
- `..` = parent directory  

- `cd <path>` = change directory  
  -- `cd ~` = go to home  
  -- `cd ..` = go up one level  

- `pwd` = print working directory  

- `ls` = list contents  
  -- `ls -a` = include hidden files  
  -- `ls -lh` = long format, human readable sizes  
  -- `ls -ltr` = list by time (oldest first)  

---

## Creating and Managing Files

- `mkdir <name>` = make directory  
- `mkdir -p path/to/dir` = make nested directories  

- `touch <file>` = create an empty file  

- `cp <source> <destination>` = copy file  
- `cp -r <dir1> <dir2>` = copy directory recursively  

- `mv <source> <destination>` = move (or rename) file  

- `nano <file>` = simple text editor inside terminal  
- `vi <file>` = more advanced text editor  

---

## Deleting Files

⚠️ **Deleted files are gone. No trash.**

- `rm <file>` = remove file  
- `rm -r <dir>` = remove directory recursively  
- `rm -rf <dir>` = force remove without asking  

---

## Viewing and Checking Files

- `cat <file>` = print whole file  
- `less <file>` = view file, scroll with arrows, quit with `q`  
- `head <file>` = first 10 lines  
- `head -n 50 <file>` = first 50 lines  
- `tail <file>` = last 10 lines  
- `tail -n 50 <file>` = last 50 lines  
- `tail -f <file>` = follow live output (logs)  

- `wc -l <file>` = count lines in file  
- `grep "pattern" <file>` = search for pattern  
- `grep -c "pattern" <file>` = count matches  

---

## Compression and Archiving

- `gzip <file.fa>` = compress (replaces file with `.gz`)  
- `gzip -c <file.fa> > file.fa.gz` = compress but keep original  
- `gzip -d <file.fa.gz>` = decompress  

- `tar -czvf output.tar.gz <dir>` = create compressed archive  
- `tar -xzvf file.tar.gz` = extract archive  

- `unzip <file.zip>` = unzip  
- `zip -r output.zip <dir>` = zip directory  

---

## System Information and Checks

- `df -h` = disk usage  
- `du -sh <dir>` = size of directory  
- `free -h` = memory usage  
- `uptime` = system load and uptime  
- `top` = interactive process monitor  
- `htop` = improved process monitor (if installed)  
- `uname -a` = kernel and system info  
- `cat /etc/os-release` = OS version  

---

## Managing Processes and Jobs

- `&` = run process in background  
- `jobs -l` = list background jobs  
- `fg` = bring background job to foreground  
- `bg` = resume job in background  
- `kill <PID>` = stop process  
- `kill -9 <PID>` = force kill process  
- `ps aux | grep <name>` = search running processes  

---

## Installing Software (with sudo)

- `sudo apt update` = update package lists  
- `sudo apt upgrade` = upgrade installed packages  
- `sudo apt install <package>` = install software  
  -- `sudo apt install build-essential`  
  -- `sudo apt install git wget curl unzip`  

- `which <command>` = check install path  
- `command -v <command>` = check if command exists  

---

## Environment and Modules

- `export PATH=$PATH:/path/to/bin` = add directory to PATH  
- `echo $PATH` = show current PATH  

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

---

## Data Transfer

- `scp <file> user@server:/path/` = copy file to server  
- `scp user@server:/path/file ./` = copy file from server  
- `rsync -avh <source> <destination>` = sync directories  
- `wget <url>` = download file from web  
- `curl -O <url>` = download file with same name  

---

## Git and Version Control

- `git clone <repo>` = clone repository  
- `git pull` = update repository  
- `git status` = check repo status  

---

## Running Pipelines (nf-core, Nextflow)

- `nextflow run nf-core/rnaseq -profile docker`  
- `nextflow run nf-core/sarek -profile conda`  
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
