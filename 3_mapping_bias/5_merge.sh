#!/bin/bash

#SBATCH --nodes=4
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=4

# Cargar módulos necesarios
module load samtools/1.9
module load hisat2/2.1.0
module load miniconda/3.9

# Variables de entrada
SNP_TAB="/home/irs16/genpsych_acc_directo/BIOBANK/RNASEQ/TFM_InesSus_RNAseq/bfiles_y_vcf/HDF5_files/snp_tab.h5"
SNP_INDEX="/home/irs16/genpsych_acc_directo/BIOBANK/RNASEQ/TFM_InesSus_RNAseq/bfiles_y_vcf/HDF5_files/snp_index.h5"
HAPLOTYPE="/home/irs16/genpsych_acc_directo/BIOBANK/RNASEQ/TFM_InesSus_RNAseq/bfiles_y_vcf/HDF5_files/haplotypes.h5"
OUTPUT_DIR="/home/irs16/genpsych_acc_directo/BIOBANK/RNASEQ/TFM_InesSus_RNAseq/bams_3_wasp"

# Archivo que contiene la lista de muestras
SAMPLES_FILE="my_samples.txt"

###################################################################################
# Paso final: Fusionar, ordenar e indexar los BAMs filtrados
for SAMPLE_NAME in $(cat ${SAMPLES_FILE}); do
    echo "Fusionando BAMs para la muestra ${SAMPLE_NAME}..."

    # Fusionar los archivos BAM filtrados
    samtools merge ${OUTPUT_DIR}/${SAMPLE_NAME}_nodup.keep.merge.bam \
                   ${OUTPUT_DIR}/${SAMPLE_NAME}_nodup.keep.bam \
                   ${OUTPUT_DIR}/${SAMPLE_NAME}_nodup.remap.keep.bam 

    echo "Ordenando el archivo fusionado ${SAMPLE_NAME}_nodup.keep.merge.bam..."

    # Ordenar el archivo BAM fusionado
    samtools sort -o ${OUTPUT_DIR}/${SAMPLE_NAME}_nodup.keep.merge.sort.bam \
                  ${OUTPUT_DIR}/${SAMPLE_NAME}_nodup.keep.merge.bam

    # Eliminar el archivo intermedio no ordenado
    rm ${OUTPUT_DIR}/${SAMPLE_NAME}_nodup.keep.merge.bam

    echo "Indexando el archivo ordenado ${SAMPLE_NAME}_nodup.keep.merge.sort.bam..."

    # Indexar el archivo BAM ordenado
    samtools index ${OUTPUT_DIR}/${SAMPLE_NAME}_nodup.keep.merge.sort.bam

    echo "Finalizado el proceso para la muestra ${SAMPLE_NAME}."
done

###################################################################################

echo "Finalizado con éxito."