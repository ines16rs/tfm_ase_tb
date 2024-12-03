#!/bin/bash

#SBATCH --nodes=4
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=4

# Cargar módulos necesarios
module load samtools/1.9
module load hisat2/2.1.0
module load miniconda/3.7

# Variables de entrada
SNP_TAB="/home/irs16/genpsych_acc_directo/BIOBANK/RNASEQ/TFM_InesSus_RNAseq/bfiles_y_vcf/HDF5_files/snp_tab.h5"
SNP_INDEX="/home/irs16/genpsych_acc_directo/BIOBANK/RNASEQ/TFM_InesSus_RNAseq/bfiles_y_vcf/HDF5_files/snp_index.h5"
HAPLOTYPE="/home/irs16/genpsych_acc_directo/BIOBANK/RNASEQ/TFM_InesSus_RNAseq/bfiles_y_vcf/HDF5_files/haplotypes.h5"
OUTPUT_DIR="/home/irs16/genpsych_acc_directo/BIOBANK/RNASEQ/TFM_InesSus_RNAseq/bams_3_wasp"

# Archivo que contiene la lista de muestras
SAMPLES_FILE="my_samples.txt"

# Filtrar lecturas remapeadas
while read SAMPLE_NAME; do
    echo "Filtrando lecturas remapeadas para la muestra ${SAMPLE_NAME}..."
    
    # Especificar los archivos de entrada y salida
    to_remap_bam=${OUTPUT_DIR}/${SAMPLE_NAME}_nodup.to.remap.bam
    remap_bam=${OUTPUT_DIR}/${SAMPLE_NAME}_nodup_remap.bam
    keep_bam=${OUTPUT_DIR}/${SAMPLE_NAME}_nodup.remap.keep.bam

    # Ejecutar el script de filtrado
    python /usr/local/wasp/0.3.4/mapping/filter_remapped_reads.py \
           ${to_remap_bam} \
           ${remap_bam} \
           ${keep_bam}

    echo "Finalizando filtrado de lecturas remapeadas para la muestra ${SAMPLE_NAME}."

done < ${SAMPLES_FILE}

###################################################################################