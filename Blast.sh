#!/bin/bash

#SBATCH -p RM-shared
#SBATCH -t 30:00:00
#SBATCH --ntasks-per-node=64
#SBATCH -J Blast

cd /ocean/projects/bio230007p/iwang1
module load BLAST
## Get Dog WGS
wget https://ftp.ncbi.nlm.nih.gov/genomes/genbank/vertebrate_mammalian/Canis_lupus/latest_assembly_versions/GCA_000002285.4_Dog10K_Boxer_Tasha/GCA_000002285.4_Dog10K_Boxer_Tasha_genomic.fna.gz
## Unzip
gzip -d GCA_000002285.4_Dog10K_Boxer_Tasha_genomic.fna.gz
## Make database using dog's database
makeblastdb -dbtype nucl -in GCA_000002285.4_Dog10K_Boxer_Tasha_genomic.fna -out dog.blastdb
## Split fasta file
awk -v size=1000 -v pre=wolf -v pad=5 '/^>/{n++;if(n%size==1){close(f);f=sprintf("%s.%0"pad"d",pre,n)}}{print>>f}' wolf.fa
## Do blast for each fa file
for i in ls wolf.*1; do
    sbatch /jet/home/iwang1/blast_single.sh $i
done
### After finished all the files, run Togff.sh