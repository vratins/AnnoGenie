#!/bin/bash
## Whole genomic sequence and RNA seq path
isPE=$1 ### Is the rna seq PE or SE
wgs=$2 ### Unzip wgs file for test species, we can make a test to check input file
rna=$3 ### Unzip rna file for related species, we can make a test to check input file
wgs_ref=$4 ## Unzip wgs file for related species
output=$5 ### Output name, string

### Test input file
if [[ "$wgs" != *".fa" ]] && [[ "$wgs" != *".fasta" ]] && [[ "$wgs" != *".fna" ]] && [[ "$wgs" != *".fq" ]] && [[ "$wgs" != *".fastq" ]]; then
        echo "Please provide fasta or fastq file for whole genomic sequence of unannotated species"
        exit
fi
if [[ "$rna" != *".fq" ]] && [[ "$rna" != *".fastq" ]]; then
        echo "Please provide fastq file for RNA Seq data of related species"
        exit
fi
if [[ "$wgs_ref" != *".fa" ]] && [[ "$wgs_ref" != *".fasta" ]] && [[ "$wgs_ref" != *".fna" ]] && [[ "$wgs_ref" != *".fq" ]] && [[ "$wgs_ref" != *".fastq" ]]; then
        echo "Please provide fasta file for whole genomic sequenceof of related species"
        exit
fi





## Build index
## Index store in the "index" folder
mkdir index 
cd index
hisat2-build -f $wgs $output
cd ../


### alignment PE
if [ $isPE=="PE" ]; then
   
else
    ## alignment SE

fi



### convert to bam and sort
samtools view -bS $output.sam > $output.bam
samtools sort $output.bam -o $output.sorted.bam

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
./stringtie/stringtie -o $output.gtf $output.sorted.bam
### Convert gtf to fa
./gffread/gffread -w $output.fa -g $wgs $output.gtf


### Make database for blast using wgs of related species
makeblastdb -dbtype nucl -in $wgs_ref -out $output.blastdb
## Blast
blastn -query $output.gtf -db $output.blastdb -out $output.temp.txt -outfmt 6 -num_threads 64
## Only select most hit
awk '!seen[$1]++'  $output.temp.txt >  $output.txt
## blast format -> gff3
if [ ! -d "genomeGTFtools" ]; then
    wget https://github.com/wrf/genomeGTFtools/archive/refs/heads/master.zip -O genomeGTFtools.zip
    unzip genomeGTFtools.zip && mv genomeGTFtools-master genomeGTFtools
fi
## To gff3
./genomeGTFtools/blast2gff.py -b $output.txt > $output.gff3