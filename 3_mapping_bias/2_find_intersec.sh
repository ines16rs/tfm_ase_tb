#!/bin/bash

#SBATCH --nodes=4
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=4

# Cargar módulos necesarios
module load python/3.7.7
module load samtools/1.9
module load hisat2/2.1.0
module load miniconda/3.7

# Variables de entrada
SNP_TAB="/home/irs16/genpsych_acc_directo/BIOBANK/RNASEQ/TFM_InesSus_RNAseq/bfiles_y_vcf/HDF5_files/snp_tab.h5"
SNP_INDEX="/home/irs16/genpsych_acc_directo/BIOBANK/RNASEQ/TFM_InesSus_RNAseq/bfiles_y_vcf/HDF5_files/snp_index.h5"
HAPLOTYPE="/home/irs16/genpsych_acc_directo/BIOBANK/RNASEQ/TFM_InesSus_RNAseq/bfiles_y_vcf/HDF5_files/haplotypes.h5"
OUTPUT_DIR="/home/irs16/genpsych_acc_directo/BIOBANK/RNASEQ/TFM_InesSus_RNAseq/bams_3_wasp"

# Crear el directorio de salida si no existe
mkdir -p ${OUTPUT_DIR}

# Archivo que contiene la lista de muestras
SAMPLES_FILE="my_samples.txt"

# Directorio donde están los BAMs
BAM_DIR="/home/irs16/genpsych_acc_directo/BIOBANK/RNASEQ/TFM_InesSus_RNAseq/bams_2_nodup"

###################################################################################
# Bloque 1:
echo "Iniciando Bloque 1"

# Iterar sobre cada muestra
while read SAMPLE_NAME; do
    echo "Identificando lecturas que pueden tener sesgos de mapeo para la muestra ${SAMPLE_NAME}..."

    # Ejecutar el script find_intersecting_snps.py para cada muestra
    python /usr/local/wasp/0.3.4/mapping/find_intersecting_snps.py \
           --is_paired_end \
           --is_sorted \
           --output_dir ${OUTPUT_DIR} \
           --snp_tab ${SNP_TAB} \
           --snp_index ${SNP_INDEX} \
           --haplotype ${HAPLOTYPE} \
           --samples ${SAMPLE_NAME} \
           ${BAM_DIR}/${SAMPLE_NAME}_nodup.bam  # El archivo BAM de entrada

    echo "Finalizando identificación de SNPs para la muestra ${SAMPLE_NAME}."

done < ${SAMPLES_FILE}
