#!/bin/bash
module load kallisto

# make index, only need to run once, input "allcds.fasta" is the fasta file including all annotated genes and unannotated ORFs.
# kallisto index -i yeast.allcdna allcdna.fasta

out=out_new/$1
mkdir -p $out
kallisto quant -i yeast.allcdna -b 100 -o $out $1_1.fastq.gz $1_2.fastq.gz 
rm $1*fastq.gz
