#!/bin/bash
#SBATCH -p RM-shared
#SBATCH -t 30:00:00
#SBATCH --ntasks-per-node=64
#SBATCH -J Blast_single

cd /ocean/projects/bio230007p/iwang1
## Merge all files
cat wolf.*1.txt > wolf_temp.txt
## Only select most hit
awk '!seen[$1]++' wolf_temp.txt > wolf.txt
rm wolf.*1*
## Get pipeline for blast format -> gff3
if [ ! -d "genomeGTFtools" ]; then
    wget https://github.com/wrf/genomeGTFtools/archive/refs/heads/master.zip -O genomeGTFtools
    unzip genomeGTFtools.zip && mv genomeGTFtools-master genomeGTFtools
fi
## blast format -> gff3
./genomeGTFtools/blast2gff.py -b wolf.txt > wolf.gff3