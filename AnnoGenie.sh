## Whole genomic sequence and RNA seq path
isPE=$1 ### Is the rna seq PE or SE
wgs=$2 ### Unzip or zip wgs file, we can make a test to check input file
rna=$3 ### Unzip or zip wgs file, we can make a test to check input file
output=$4 ### Output name, string



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

./stringtie/stringtie -o $output.gtf $output.sorted.bam
./gffread/gffread -w $output.fa -g $wgs $output.gtf
