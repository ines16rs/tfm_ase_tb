# Análisis de expresión alélica específica en trastorno bipolar

Este TFM se centra en investigar la expresión alélica específica (ASE) en individuos con trastorno bipolar (TB), una enfermedad psiquiátrica compleja, cuya base genética todavía no se conoce por completo. Como datos de partida para este análisis se han utilizado las lecturas de RNA-seq y los genotipos de individuos afectados por TB y controles sanos, para tratar de identificar factores genéticos potencialmente implicados en el desarrollo de esta enfermedad. Para ello, se han utilizado herramientas bioinformáticas como **WASP** y **AllelicImbalance**, para tratar de identificar posibles diferencias de expresión alélica. 

---

## Contenido del este repositorio

A continuación, se detallan las carpetas que forman este repositorio, junto con una breve descripción de los archivos que contiene cada una de ellas:

### 0_bfiles_to_vcf
- **0_bfiles_to_VCFmultisample:** conversión de archivos bfiles a formato VCF (multimuesta).
- **1_generar_22vcfs:** generación de un archivo VCF para cada cromosoma autosómico (1 al 22) a partir del VCF multimuestra completo.

### 1_generacion_bams
- **0_fastqc:** control de calidad de las lecturas de RNA-seq empleando **FastQC**.
- **1_alignment:** alineamiento de las lecturas de RNA-seq al genoma de referencia (hg38) utilizando la herramienta **hisat2**.
- **Run1_alignment:** script para ejecutar el alineamiento en el centro de computación científica de la UAM.

### 2_markdup_samtools
- **2_markdup:** eliminación aleatoria de duplicados utilizando **samtools markdup**.

### 3_mapping_bias
- **0_generar_my_sample:** generación de un archivo de texto con los identificadores de las muestras incluidas en el estudio.
- **1_VCF_to_HDF5:** creación de los archivos **HDF5** necesarios para **WASP** a partir de los VCFs de los 22 cromosomas.
- **2_find_intersec:** identificación de las lecturas con sesgos de mapeo potenciales usando **find_intersecting_snps.py** de WASP.
- **3_remapeo:** remapeo de lecturas con sesgos potenciales utilizando **hisat2**.
- **4_filtrar_remap:** filtrado de lecturas remapeadas con **filter_remapped_reads.py** de WASP.
- **5_merge:** unión, ordenación e indexación de las lecturas que han superado los filtros anteriores, utilizando **samtools**.

### 4_allelicimbalanceR
- **1_allelicimbalance_binom:** análisis de ASE con **AllelicImbalance** en R: instalación de paquetes, recuento alélico, filtrado de variantes heterocigotas y lecturas mínimas, y test binomial por variante e individuo.

### 5_wilcoxon_rank_sum
- **1_proporciones_y_wilcoxon:** cálculo de las proporciones del alelo de referencia y aplicación del test de Wilcoxon Rank Sum para comparar la ASE de los dos grupos.

### 6_fisher_exact
- **1_tablas_conting_y_fisher:** generación de las tablas de contingencia y aplicación del test exacto de Fisher para analizar los patrones de expresión mayoritaria.
- **2_generar_tabla_resumen_info:** creación de una tabla resumen con información clave de cada variante analizada (rsID, posición, p-valores, genes, etc.).

### 7_graficos_y_tablas
- **1_generacion_boxplots_bams:** creación de boxplots del número de lecturas RNA-seq en las diferentes etapas de procesamiento, separadas por condición clínica.
- **2_tabla_n_variantes_filtros:** código para generar una tabla con la distribución de variantes en cada etapa del filtrado.
- **3_manhattan_plot_variantes:** generación de un gráfico de tipo **Manhattan Plot** para las variantes analizadas en el test de Wilcoxon Rank Sum.
- **4_locuszoom_link:** enlace de acceso a la herramienta **LocusZoom** para visualizar los datos de los genes comunes entre el análisis de ASE y el de *"gene-based"*.

### 8_anotacion_y_enriq_funcional
- **0_external_gene_IDs_enrichR:** extracción de **external_gene_IDs** para el análisis de enriquecimiento funcional con **Enrichr**.
- **1_Enrichr_link:** enlace a la herramienta **Enrichr** para el análisis de enriquecimiento funcional.
- **2_VEP_Ensembl_link:** enlace a **VEP** de Ensembl para anotación funcional (valores **CADD** y frecuencias alélicas de **gnomAD**).

### Manuales_y_documentación
- Subcarpetas con la documentación recopilada de las principales herramientas utilizadas en el desarrollo del TFM.

### Viñeta TFM
- Viñeta explicativa del trabajo desarrollado en este TFM (pdf y html)
---

