#!/bin/bash

module load samtools/1.16.1

# Directorios de entrada y salida
input_dir="/home/irs16/genpsych_acc_directo/BIOBANK/RNASEQ/TFM_InesSus_RNAseq/bams_1_iniciales"
output_dir="/home/irs16/genpsych_acc_directo/BIOBANK/RNASEQ/TFM_InesSus_RNAseq/bams_2_nodup"

# Crear el directorio de salida si no lo hemos creado previamente
mkdir -p "$output_dir"

# Iterar sobre cada archivo .bam en el directorio de entrada
for bam_file in "$input_dir"/*.bam; do
    # Obtener el nombre base del archivo BAM sin la extensión .bam
    base_name=$(basename "$bam_file" .bam)

    # Definir los nombres de los archivos intermedios y de salida
    queryname_sorted_bam="$output_dir/${base_name}_queryname_sorted.bam"
    fixmate_bam="$output_dir/${base_name}_fixmate.bam"
    sorted_bam="$output_dir/${base_name}_sorted.bam"
    final_bam="$output_dir/${base_name}_nodup.bam"

    # Paso 1: Ordenar por nombre de consulta (queryname)
    echo "Ordenando por nombre de consulta $bam_file"
    samtools sort -n -o "$queryname_sorted_bam" "$bam_file" 2>> "$output_dir/sort_queryname_errors.log"

    # Paso 2: Ejecutar samtools fixmate y registrar errores
    echo "Ejecutando samtools fixmate para $queryname_sorted_bam"
    samtools fixmate -m "$queryname_sorted_bam" "$fixmate_bam" 2>> "$output_dir/fixmate_errors.log"

    # Paso 3: Ejecutar samtools sort por coordenadas y registrar errores
    echo "Ejecutando samtools sort por coordenadas para $fixmate_bam"
    samtools sort -o "$sorted_bam" "$fixmate_bam" 2>> "$output_dir/sort_errors.log"

    # Paso 4: Ejecutar samtools markdup (con -r para "remove" aleatoriamente) y registrar errores
    echo "Ejecutando samtools markdup para $sorted_bam"
    samtools markdup -r "$sorted_bam" "$final_bam" 2>> "$output_dir/markdup_errors.log"

    # Paso 5: Generar el índice .bam.bai para el nuevo BAM sin duplicados
    echo "Generando índice para $final_bam"
    samtools index "$final_bam"

    # Paso 6: Generar estadísticas rápidas con flagstat para el BAM inicial
    echo "Generando estadísticas iniciales con flagstat para $bam_file"
    samtools flagstat "$bam_file" > "$output_dir/${base_name}_initial_flagstat.txt"

    # Paso 7: Generar estadísticas rápidas con flagstat para el BAM sin duplicados
    echo "Generando estadísticas con flagstat para $final_bam"
    samtools flagstat "$final_bam" > "$output_dir/${base_name}_nodup_flagstat.txt"

    # Limpiar archivos intermedios si no son necesarios
    rm "$queryname_sorted_bam"
    rm "$fixmate_bam"
    rm "$sorted_bam"

done

echo "Proceso completado. Se pueden revisar los logs de errores y las estadísticas en la carpeta $output_dir"
