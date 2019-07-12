#!/bin/bash
module load kallisto
out=out_new/$1
mkdir -p $out
[[ -f $1_1.fastq.gz && -f $1_2.fastq.gz ]] && kallisto quant -i yeast.allcds -b 100 -o $out $1_1.fastq.gz $1_2.fastq.gz || kallisto quant -i yeast.allcds --single -l 200 -s 20 -b 100 -o $out $1*.fastq.gz
rm $1*fastq.gz
