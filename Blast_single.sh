#!/bin/bash
#SBATCH -p RM-shared
#SBATCH -t 30:00:00
#SBATCH --ntasks-per-node=64
#SBATCH -J Blast_single

cd /ocean/projects/bio230007p/iwang1
module load BLAST
file=$1
echo "start blast $file"
blastn -query $file -db dog.blastdb -out $file.txt -outfmt 6



