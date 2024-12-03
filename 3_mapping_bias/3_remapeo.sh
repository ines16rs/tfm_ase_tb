#!/bin/bash

# Cargar módulos necesarios
module load samtools/1.9
module load hisat2/2.1.0

# Variables de entrada
SNP_TAB="/home/irs16/genpsych_acc_directo/BIOBANK/RNASEQ/TFM_InesSus_RNAseq/bfiles_y_vcf/HDF5_files/snp_tab.h5"
SNP_INDEX="/home/irs16/genpsych_acc_directo/BIOBANK/RNASEQ/TFM_InesSus_RNAseq/bfiles_y_vcf/HDF5_files/snp_index.h5"
HAPLOTYPE="/home/irs16/genpsych_acc_directo/BIOBANK/RNASEQ/TFM_InesSus_RNAseq/bfiles_y_vcf/HDF5_files/haplotypes.h5"
OUTPUT_DIR="/home/irs16/genpsych_acc_directo/BIOBANK/RNASEQ/TFM_InesSus_RNAseq/bams_3_wasp"

# Remapeo de lecturas con potencial sesgo de mapeo
for remap_fq1 in ${OUTPUT_DIR}/*.remap.fq1.gz; do
    # Extraer el nombre de la muestra del nombre del archivo
    SAMPLE_NAME=$(basename ${remap_fq1} .remap.fq1.gz)
    echo "Nombre de la muestra ${SAMPLE_NAME} extraído con éxito."

    # Especificar el archivo de lectura emparejada
    remap_fq2=${OUTPUT_DIR}/${SAMPLE_NAME}.remap.fq2.gz

    echo "Remapeando lecturas para la muestra ${SAMPLE_NAME}..."

    # Ejecutar hisat2 para remapear las lecturas emparejadas
    hisat2 --rna-strandness RF -k 1 -p 4 -x /home/irs16/genpsych_acc_directo/BIOBANK/RNASEQ/TFM_InesSus_RNAseq/fasta_index/hg38 \
           -1 ${remap_fq1} \
           -2 ${remap_fq2} \
           -S ${OUTPUT_DIR}/${SAMPLE_NAME}_remap.sam 2>> summary_remapping_${SAMPLE_NAME}.txt

    echo "Convirtiendo ${SAMPLE_NAME}_remap.sam a BAM..."

    # Convertir a BAM
    samtools view -@ 4 -bo ${OUTPUT_DIR}/${SAMPLE_NAME}_remap_notsorted.bam ${OUTPUT_DIR}/${SAMPLE_NAME}_remap.sam

    echo "Eliminando ${SAMPLE_NAME}_remap.sam"
    rm ${OUTPUT_DIR}/${SAMPLE_NAME}_remap.sam

    echo "Ordenando ${SAMPLE_NAME}_remap_notsorted.bam..."
    samtools sort -@ 4 -o ${OUTPUT_DIR}/${SAMPLE_NAME}_remap.bam ${OUTPUT_DIR}/${SAMPLE_NAME}_remap_notsorted.bam

    echo "Indexando ${SAMPLE_NAME}_remap.bam..."
    samtools index ${OUTPUT_DIR}/${SAMPLE_NAME}_remap.bam

    echo "Finalizando remapeo de lecturas para la muestra ${SAMPLE_NAME}."
done

###################################################################################