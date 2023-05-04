#!/bin/bash
## Whole genomic sequence and RNA seq path
wkdir=$1 ### name of the working directory, string
isPE=$2 ### Is the rna seq PE or SE
wgs=$3 ### path for unzipped wgs file for test species, we can make a test to check input file
rna=$4 ### path for directory storing unzipped rna file for related species, we can make a test to check input file
wgs_ref=$5 ## path for unzipped wgs file for related species
related=$6 ## name of the related species
output=$7 ### Output name for unannotated species, string

### sample command for downloading pair-end rna seq file
### fastq-dump --outdir [output directory] --split-files [ID]



### Test input file
if [[ "$wgs" != *".fa" ]] && [[ "$wgs" != *".fasta" ]] && [[ "$wgs" != *".fna" ]] && [[ "$wgs" != *".fq" ]] && [[ "$wgs" != *".fastq" ]]; then
        echo "Please provide fasta or fastq file for whole genomic sequence of unannotated species"
        exit
fi

if [[ "$wgs_ref" != *".fa" ]] && [[ "$wgs_ref" != *".fasta" ]] && [[ "$wgs_ref" != *".fna" ]] && [[ "$wgs_ref" != *".fq" ]] && [[ "$wgs_ref" != *".fastq" ]]; then
        echo "Please provide fasta file for whole genomic sequenceof of related species"
        exit
fi

# Check if the rna argument is a directory
if [[ ! -d "$rna" ]]; then
  echo "$rna is not a directory, plase provide the directory storing the RNAseq files"
  exit
fi

# Count the number of files in the rna directory
count=$(ls -1 "$rna" | wc -l)

# Check if input file and isPE option is compatible
if [[ $count -eq 1 ]] && [[ $isPE=="PE" ]]; then
  echo "The directory $rna contains one file but the mode is Pair-end."
  exit
elif [[ $count -eq 2 ]] && [[ $isPE=="SE" ]]; then
  echo "The directory $rna contains two files but the mode is Pair-end."
  exit
else
  echo "The directory $rna contains $count files, proceed to next test."
fi


# Check if all files in the directory end with ".fq" or ".fastq"
if [[ $(find "$rna" -maxdepth 1 -type f ! -name "*.fq" ! -name "*.fastq" | wc -l) -ne 0 ]]; then
  echo "Not all files in $rna end with .fq or .fastq"
  exit 1
fi

echo "All files in $rna end with .fq or .fastq, proceed to fastQC"



### make working directory
mkdir $wkdir
cd $wkdir


# Get the list of file names that end with ".fq" or ".fastq"
fq_files=$(find "$1" -maxdepth 1 -type f \( -name "*.fq" -o -name "*.fastq" \))

if [ $isPE=="PE" ]; then
    file1=$(echo "$fq_files" | head -n 1)
    file2=$(echo "$fq_files" | head -n 2 | tail -n 1)
else 
    file1=$(echo "$fq_files" | head -n 1)

### fastqc

mkdir QC
cd QC

if [ $isPE=="PE" ]; then
    fastqc $file1 $file2
else 
    fastqc $file1 

cd ..

mkdir trim
cd trim
### trimmomatic
if [ $isPE=="PE" ]; then
    trimmomatic PE $file1 $file2 output_forward_paired.fastq.gz output_forward_unpaired.fastq.gz output_reverse_paired.fastq.gz output_reverse_unpaired.fastq.gz ILLUMINACLIP:NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    gzip -d output_forward_paired.fastq.gz 
    gzip -d output_reverse_paired.fastq.gz
else 
    trimmomatic SE $file1 output.fastq.gz ILLUMINACLIP:NexteraPE-PE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:32
    gzip -d output.fastq.gz 


cd ..

### fastqc check after trim

mkdir QC2
cd QC2

if [ $isPE=="PE" ]; then
    fastqc ../trim/output_forward_paired.fastq ../trim/output_reverse_paired.fastq
else 
    fastqc ../trim/output.fastq

cd ..

mkdir align
cd align
## Build index
## Index store in the "index" folder
mkdir index 
cd index
hisat2-build -f $wgs $related
cd ..


if [ $isPE=="PE" ]; then
    ### alignment PE
    hisat2 -q -x ./index/$related -1 ../trim/output_forward_paired.fastq -2 ../trim/output_reverse_paired.fastq -S $output.sam

else
    ### alignment SE
    hisat2 -q -x ./index/$related -U ../trim/output.fastq -S $output.sam

fi


### convert to bam and sort
samtools view -bS $output.sam > $output.bam
samtools sort $output.bam -o $output.sorted.bam

cd ..

### Clone stringtie
if [ ! -d "stringtie" ]; then
    git clone https://github.com/gpertea/stringtie
    cd stringtie
    make release
    cd ../
    echo "StringTie installed"
fi

### Clone gffread
if [ ! -d "gffread" ]; then
    git clone https://github.com/gpertea/gffread
    cd gffread
    make release
    cd ../
    echo "gffread installed"
fi

### Transcripts
./stringtie/stringtie -o $output.gtf ./align/$output.sorted.bam
### Convert gtf to fa
./gffread/gffread -w $output.fa -g $wgs $output.gtf


### Make database for blast using wgs of related species
makeblastdb -dbtype nucl -in $wgs_ref -out $output.blastdb
## Blast
blastn -query $output.fa -db $output.blastdb -out $output.temp.txt -outfmt 6 -num_threads 64
## Only select most hit
awk '!seen[$1]++'  $output.temp.txt >  $output.txt
## blast format -> gff3
if [ ! -d "genomeGTFtools" ]; then
    wget https://github.com/wrf/genomeGTFtools/archive/refs/heads/master.zip -O genomeGTFtools.zip
    unzip genomeGTFtools.zip && mv genomeGTFtools-master genomeGTFtools
fi
## To gff3
./genomeGTFtools/blast2gff.py -b $output.txt > $output.gff3