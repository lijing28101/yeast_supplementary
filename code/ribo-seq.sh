#!/bin/bash

#./ribo-seq.sh SRRXXXX

#input directory containing SRR files

# Load all required modules, files, prereqs. These are
# Modules: HISAT2, aspera, samtools
# ribotricer, bbmap, and sratoolkit 2.9.6 already download
# Files: yeastgenome.fa, corresponding yeastAnnotation.gtf, SRR files, rRNA seq 

module load hisat2
module load aspera
module load samtools

openssh=/shared/hpc/aspera/cli/3.7.7/etc/asperaweb_id_dsa.openssh
bbmap_folder=/home/jingli/bbmap/
bbduk_ref=${bbmap_folder}/resources/adapters.fa
hisat2_index=/ptmp/LAS/jingli/yeast_ribo/hisat2_index/histat2_index
rRNA_ref=/ptmp/LAS/jingli/yeast_ribo/rRNA.fa
bam_out=/ptmp/LAS/jingli/yeast_ribo/bam_out
ribotricer_index=/ptmp/LAS/jingli/yeast_ribo/yeast_riboindex_candidate_orfs.tsv
ribo_out=/ptmp/LAS/jingli/yeast_ribo/ribo_out

#input srr sample ID.
line=$1

#Download SRA data
srr=$(echo "${line::3}" | tr '[:upper:]' '[:lower:]')
num=$(echo "${line: -1}")
part=$(echo "${line::6}")

if [[ $(echo -n $line | wc -m) -eq 10 ]]; then
  file=$(echo "era-fasp@fasp.sra.ebi.ac.uk:/vol1/${srr}/${part}/00${num}/${line}")
else
  file=$(echo "era-fasp@fasp.sra.ebi.ac.uk:/vol1/${srr}/${part}/${line}")
fi

ascp -i ${openssh} -P33001 -QT -l 500m ${file} ${line}.sra

if [ $? -ne 0 ]; then
            echo "FAILED TO DOWNLOAD " "${line}"
else
            echo "Download successful " "${line}"
fi

#Convert sra to fastq file
/work/LAS/mash-lab/jing/bin/sratoolkit.2.9.6-1-centos_linux64/bin/fasterq-dump  -f -o ${line}_pass_1.fastq ${line}.sra


#delete adapter sequence
${bbmap_folder}/bbduk.sh in=${line}.fastq out=${line}_clean.fq ref=${bbduk_ref} ktrim=r k=11 mink=4 hdist=1 qtrim=rl trimq=10
#Some old samples may use sanger sequencing, then add qin=33
#${bbmap_folder}/bbduk.sh in=${line}.fastq out=${line}_clean.fq ref=${bbduk_ref} qin=33 ktrim=r k=11 mink=4 hdist=1 qtrim=rl trimq=10

#delete rRNA reads
${bbmap_folder}/bbsplit.sh in=${line}_clean.fq ref=${rRNA_ref} basename=out_${line}_%.fq outu=${line}_bbspout.fq


#alignment by hisat
#build hisat2 index 
#hisat2-build Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa hisat2_index/histat2_index

hisat2 --mp 1,1 --no-spliced-alignment --end-to-end --rdg 10000,10000 --rfg 10000,10000 -p 28 -x ${hisat2_index} -U ${line}_bbspout.fq -S ${line}_hisat.sam

#convert sam to bam, and delete old sam file
samtools view -@ 28 -b -o ${line}_hisat.bam ${line}_hisat.sam
samtools sort -o ${bam_out}/${line}_sorted.bam  -T SRR6234838_temp --threads 28 ${line}_hisat.bam
rm ${line}_hisat.sam ${line}_hisat.bam

#find translating genes and ORFs
#build ribotricer index
# ribotricer prepare-orfs --gtf Saccharomyces_cerevisiae.updated.gtf --fasta Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa --prefix yeast_riboindex --min_orf_length 18
ribotricer detect-orfs --bam ${bam_out}/${line}_sorted.bam --ribotricer_index ${ribotricer_index} --prefix ${ribo_out}/${line} --phase_score_cutoff 0.3 --min_valid_codons 3
