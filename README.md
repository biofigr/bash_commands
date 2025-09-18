# bash_commands

Create a list of common bash commands for biologists new to CLI in one place.

## Navigating CLI

- / = root
- ~ = home
- . = current directory
- .. = previous directory
- cd = change directory

- pwd = print working directory
- ls = list
- ls -a = list all (including hidden files)
- ls -lh = list long, human readable

## Creating files

- mkdir = make directory
- cp = copy
- mv = move

## Deleting files

### N.B. There is no trash folder, careful what you delete!

- rm = remove file
- rm -rf = remove folder and contents

## Compression

- gzip -c <file.fa> > <file.fa.gz> = compress, but keep original file untouched
- gzip -d = decompress
- tar -czvf <output.tar.gz> /path/to/folder/ = compress (c) with gzip (z), verbose as it runs (v) and file name (f)
- tar -xzvf file.tar.gz = extract (x)

## System checks

- df -h = show disk space
- free -h = show memory
- top = show jobs

## Kill a job

- jobs -l = list jobs
# [1]+ 37711 Running                 nohup tar -czf all_fastq.tar.gz ./data/ &
kill 37711
