#!/bin/bash

# Define el path donde están las carpetas individuales con los archivos .fq.gz
# data_path="/home/irs16/genpsych_acc_directo/BIOBANK/RNASEQ/TFM_InesSus_RNAseq/fastqs"

# Cambiar al directorio donde están las carpetas individuales
# cd $data_path

# Ejecutar el script de alineación para cada carpeta
for carpeta in */; do
    # Obtener el nombre de la carpeta (sin la barra diagonal final)
    sample_name=${carpeta%/}

    # Buscar el archivo forward que empiece por el nombre de la muestra y termine en _1.fq.gz
    forward=$(find "$carpeta" -name "${sample_name}*_1.fq.gz" -print -quit)
    # Buscar el archivo reverse que empiece por el nombre de la muestra y termine en _2.fq.gz
    reverse=$(find "$carpeta" -name "${sample_name}*_2.fq.gz" -print -quit)

    # Verificar si los archivos forward y reverse existen
    if [[ -f "$forward" && -f "$reverse" ]]; then
        # Enviar el trabajo a la cola con los argumentos
        sbatch -A genpsych_serv -p bioinfo 1_alignment.sh "$sample_name" "$forward" "$reverse"
    else
        echo "Archivos no encontrados para la muestra: $sample_name"
    fi
done

