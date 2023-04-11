#!/bin/bash
#SBATCH -p RM-shared
#SBATCH -t 03:00:00
#SBATCH --ntasks-per-node=64
#SBATCH -J alignment

module load hisat2
module load samtools/1.13.0

cd /ocean/projects/bio230007p/lqiao1/alignment

### alignment
hisat2 -x /ocean/projects/bio230007p/lqiao1/alignment/index/wolf -1 /ocean/projects/bio230007p/lqiao1/RNA-seq/sample/ERR10821899_1.fastq -2 /ocean/projects/bio230007p/lqiao1/RNA-seq/sample/ERR10821899_2.fastq -S ERR10821899.sam

### convert to bam and sort
samtools view -bS ERR10821899.sam > ERR10821899.bam

samtools sort ERR10821899.bam -o ERR10821899.sorted.bam