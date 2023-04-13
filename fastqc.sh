module load sra-tools
##Download the RNA-seq files from the SRA id into specified directory (Paired-end sequences)
fastq-dump --outdir /ocean/projects/bio230007p/vshende/rnaseq --split-files DRR453607
##Run fastqc to check quality of data
#fastqc file1 file2

trimmomatic PE DRR453607_1.fastq DRR453607_2.fastq -baseout output.fq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36