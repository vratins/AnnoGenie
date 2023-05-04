#!/bin/bash

awk 'BEGIN{FS=OFS="\t"} {
    # extract fields
    split($1, a, "_"); chr=a[5]
    split($1, a, "_"); gene = a[1] "_" a[2] "_" a[3]
    start=$4
    end=$5
    score=$6
    plus = $7
    sub(/GeneID:/, "", gene)
    sub(/\..*$/, "", gene)
    
    # convert to 1-based half-open format
    if ($4 == "+") {
        start += 1
    } else {
        end -= 1
    }

    # # convert chromosome name to UCSC format (if needed)
    # if (chr ~ /^chr/) {
    #     sub(/^chr/, "", chr)
    # } else {
    #     chr = "chr" chr
    # }
    
    # output modified line
    # print chr"\tblastn\tgene\t"start"\t"end"\t"\t"$4"\t\tgene_id="gene
    # print "$chr\tblastn\tgene\t$start\t$end\t.\t$4\t.\tgene_id=$gene"
    OFS="\t"; print chr, ".", "gene", start, end, score, plus, ".", "gene_id=" gene
}' ./wolf_1.gff3 > output3.gff3
