# Análisis de expresión alélica específica en trastorno bipolar

Este TFM se centra en investigar la expresión alélica específica (ASE) en individuos con trastorno bipolar (TB), una enfermedad psiquiátrica compleja cuya base genética aún no se conoce por completo. Se han utilizado datos de RNA-seq y genotipos de individuos afectados por TB y controles sanos, con el objetivo de identificar factores genéticos implicados en la etiología y patogénesis de esta enfermedad. Este análisis podría contribuir a mejorar el diagnóstico y tratamiento del TB. Para llevar a cabo el estudio, se han empleado herramientas bioinformáticas como **WASP** y **AllelicImbalance**.

---

## Estructura del Repositorio

### 0_bfiles_to_vcf
- **0_bfiles_to_VCFmultisample:** Conversión de archivos bfiles a formato VCF multimuesta.
- **1_generar_22vcfs:** Generación de un archivo VCF para cada cromosoma autosómico a partir del VCF total.

### 1_generacion_bams
- **0_fastqc:** Control de calidad de las lecturas de RNA-seq.
- **1_alignment:** Alineamiento de las lecturas de RNA-seq al genoma de referencia (hg38).
- **Run1_alignment:** Script para ejecutar el alineamiento en el centro de computación científica de la UAM.

### 2_markdup_samtools
- **2_markdup:** Eliminación aleatoria de duplicados utilizando **samtools markdup**.

### 3_mapping_bias
- **0_generar_my_sample:** Generación de un archivo de texto con los identificadores de las muestras del estudio.
- **1_VCF_to_HDF5:** Creación de archivos **HDF5** necesarios para **WASP** a partir de los VCFs de los 22 cromosomas.
- **2_find_intersec:** Identificación de lecturas con sesgos de mapeo usando **find_intersecting_snps.py** de WASP.
- **3_remapeo:** Remapeo de lecturas con sesgos utilizando **hisat2**.
- **4_filtrar_remap:** Filtrado de lecturas remapeadas con **filter_remapped_reads.py** de WASP.
- **5_merge:** Fusión, ordenación e indexación de los BAMs filtrados utilizando **samtools**.

### 4_allelicimbalanceR
- **1_allelicimbalance_binom:** Análisis de ASE con **AllelicImbalance** en R: instalación de paquetes, recuento alélico, filtrado y test binomial por variante e individuo.

### 5_wilcoxon_rank_sum
- **1_proporciones_y_wilcoxon:** Cálculo de proporciones del alelo de referencia y aplicación del test de Wilcoxon Rank Sum para comparar los grupos.

### 6_fisher_exact
- **1_tablas_conting_y_fisher:** Generación de tablas de contingencia y test exacto de Fisher para patrones de expresión.
- **2_generar_tabla_resumen_info:** Creación de una tabla resumen con información clave de cada variante (rsID, posición, p-valores, genes, etc.).

### 7_graficos_y_tablas
- **1_generacion_boxplots_bams:** Creación de boxplots del número de lecturas RNA-seq en diferentes etapas de procesamiento, separadas por condición clínica.
- **2_tabla_n_variantes_filtros:** Tabla con la distribución de variantes en cada etapa del filtrado.
- **3_manhattan_plot_variantes:** Generación de un **Manhattan Plot** para las variantes analizadas en el test de Wilcoxon.
- **4_locuszoom_link:** Acceso a **LocusZoom** para visualizar genes comunes en diferentes análisis.

### 8_anotacion_y_enriq_funcional
- **0_external_gene_IDs_enrichR:** Extracción de **external_gene_IDs** para el análisis de enriquecimiento funcional con **Enrichr**.
- **1_Enrichr_link:** Enlace a la herramienta **Enrichr** para el análisis de enriquecimiento funcional.
- **2_VEP_Ensembl_link:** Enlace a **VEP** de Ensembl para anotación funcional (valores **CADD** y frecuencias **gnomAD**).

---

