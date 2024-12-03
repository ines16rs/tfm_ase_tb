#!/bin/bash
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=4
module load hisat2/2.1.0
module load samtools/1.9
index=$1
forward=$2
reverse=$3

#Align each sample
echo Starting with sample ${index}
hisat2 --rna-strandness RF -k 1 -p 4 -x /home/irs16/genpsych_acc_directo/BIOBANK/RNASEQ/TFM_InesSus_RNAseq/fasta_index/hg38 -1 ${forward} -2 ${reverse} -S ${index}.sam 2>> summary_alignment_${index}.txt

echo Getting into samtools, sample ${index}
samtools view -@ 4 -bo ${index}_notsorted.bam ${index}.sam

echo Removing ${index}.sam
rm ${index}.sam

echo Sorting ${index}_notsorted.bam
samtools sort -@ 4 -o ${index}.bam ${index}_notsorted.bam

echo Indexing
samtools index ${index}.bam

echo Finishing with sample ${index}


