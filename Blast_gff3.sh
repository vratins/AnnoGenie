#!/bin/bash

#SBATCH -p RM-shared
#SBATCH -t 30:00:00
#SBATCH --ntasks-per-node=64
#SBATCH -J Blast_All

### Script for runing Blast without splitting file, and converts output to gff3
## For dog/wolf takes 3~4 hr
cd /ocean/projects/bio230007p/iwang1
module load BLAST
## Get Dog WGS
wget https://ftp.ncbi.nlm.nih.gov/genomes/genbank/vertebrate_mammalian/Canis_lupus/latest_assembly_versions/GCA_000002285.4_Dog10K_Boxer_Tasha/GCA_000002285.4_Dog10K_Boxer_Tasha_genomic.fna.gz
## Unzip
gzip -d GCA_000002285.4_Dog10K_Boxer_Tasha_genomic.fna.gz
## Make database using dog WGS
makeblastdb -dbtype nucl -in GCA_000002285.4_Dog10K_Boxer_Tasha_genomic.fna -out dog.blastdb
## Blast
blastn -query wolf.fa -db dog.blastdb -out wolf_temp.txt -outfmt 6 -num_threads 64
## Only select most hit
awk '!seen[$1]++' wolf_temp.txt > wolf.txt
## Get pipeline for blast format -> gff3
if [ ! -d "genomeGTFtools" ]; then
    wget https://github.com/wrf/genomeGTFtools/archive/refs/heads/master.zip -O genomeGTFtools.zip
    unzip genomeGTFtools.zip && mv genomeGTFtools-master genomeGTFtools
fi
## To gff3
./genomeGTFtools/blast2gff.py -b wolf.txt > wolf.gff3