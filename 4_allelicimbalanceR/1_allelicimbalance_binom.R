setwd("C:/irs_tfm")
getwd()

################################################################################
### 1. Instalar los paquetes necesarios (AllelicImbalance y VariantAnnotation)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("AllelicImbalance", "VariantAnnotation",
                       "GenomicAlignments", "SNPlocs.Hsapiens.dbSNP144.GRCh38"))

install.packages("readxl")

if (!requireNamespace("biomaRt", quietly = TRUE)) {
  install.packages("biomaRt")
}

# Instalar dplyr si aún no lo tienes
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}

# Cargar AllelicImbalance, VariantAnnotation y otros paquetes necesarios
library(AllelicImbalance)
library(VariantAnnotation)
library(GenomicAlignments)
library(readxl)
library(SNPlocs.Hsapiens.dbSNP144.GRCh38)
library(biomaRt)
library(dplyr)
library(tidyr)

################################################################################
### 2. Rutas a los archivos necesarios:
# Definir la ruta al archivo VCF (VCF multisample)
pathToVcf <- "C:/irs_tfm/vcf_multi.vcf"

# Definir el directorio donde se encuentran los BAMs 
pathToFiles <- "X:/ASE_irs/bams_ase"

# Se carga el excel con los fenotipos de los individuos (affected vs unaffected)
fenotipo <- read_excel("C:/irs_tfm/fenotipo_cond.xlsx")

################################################################################
### 3. Leer el archivo VCF (VCF: info de todos los individuos) y lo filtramos
VCF_raw <- readVcf(pathToVcf, "hg38") 
# Eliminamos los posibles duplicados y el archivo con duplicados para liberar espacio
VCF <- unique(VCF_raw)
rm(VCF_raw)
cat("Archivo VCF cargado exitosamente.\n")

# Convertir el VCF en un objeto GRanges
gr <- granges(VCF)

# Filtrado: quedarnos solo con las variantes bialélicas (Alelos: 1 REF y 1 ALT)
gr.filt <- gr[width(mcols(gr)[,"REF"]) == 1 & 
                sapply(mcols(gr)[,"ALT"], function(alt) length(alt) == 1 &&
                         !is.na(alt) && nchar(alt) > 0)]
save(gr.filt, file = "gr.filt.RData")

cat("Filtrado completado: se han retenido solo variantes bialélicas.\n")

# Filtrar directamente los primeros 22 niveles de gr.filt_chr (chr 1 al 22)
gr.filt <- keepSeqlevels(gr.filt, seqlevels(gr.filt)[1:22], pruning.mode = "coarse")

# Filtrar variantes bialélicas: REF de longitud 1 y ALT con exactamente un alelo
# Comprobar si todas las variantes son bialélicas (en todas las posiciones)
biallelic_check <- sapply(mcols(gr.filt)$ALT, function(alt) length(alt) == 1)
result_bialelic_check <- sum(biallelic_check) == length(gr.filt) 
cat("Verificación de variantes bialélicas completada: ", result_bialelic_check, "\n")
# Debería salir TRUE. Sí, sale TRUE

# Guardar el resultado en un txt
save(result_bialelic_check, file = "result_bialelic_check.Rdata")

# Función readGT: permite leer la info de fase del archivo VCF de partida
geno_fullVCF <- geno(VCF)$GT
save(geno_fullVCF, file = "geno_fullVCF.Rdata")
cat("Información genotípica del VCF leida correctamente. \n")

################################################################################
### 4. Trabajo con un solo cromosoma en cada iteración:

# Definir los tamaños de los cromosomas humanos (autosómicos) (versión hg38)
chr_sizes <- c(248956422, 242193529, 198295559, 190214555, 181538259, 
               170805979, 159345973, 145138636, 138394717, 133797422, 
               135086622, 133275309, 114364328, 107043718, 101991189, 
               90338345, 83257441, 80373285, 58617616, 64444167, 
               46709983, 50818468) #Se sacan de SNPlocs.Hsapiens.dbSNP144.GRCh38

# Crear un bucle que itere a través de los cromosomas 1 al 22 (autosómicos)
for (chr in 5:1) {
  
  # Definir el nombre del cromosoma actual y la longitud del cromosoma
  chr_name <- paste0("chr", chr)
  chr_length <- chr_sizes[chr]
  
  # Crear un directorio específico para guardar los resultados del chr actual
  result_dir <- paste0("resultados_", chr_name) # Directorio para guardar resultados
  dir.create(result_dir, showWarnings = FALSE)
  
  # Informar de que el proceso ha comenzado para el cromosoma actual
  cat(paste0("Iniciando análisis para ", chr_name, "...\n"))
  
  ############################################
  ### 4.1. Filtrar variantes del cromosoma estudiado del gr o vcf total
  gr.filt_chr <- gr.filt[seqnames(gr.filt) == c(chr)]
  
  # Eliminamos los posibles duplicados
  gr.filt_chr <- gr.filt_chr[!duplicated(gr.filt_chr)]
  
  ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######
  
  # Filtrar para conservar solo los cromosomas 1 al 22 en el objeto gr_filt_chr
  chromosomes_to_keep <- paste0("chr", 1:22)  # Lista de los cromosomas del 1 al 22
  
  # Crear un nuevo objeto con solo estos cromosomas
  gr.filt_chr <- keepSeqlevels(gr.filt_chr, 1:22, pruning.mode = "coarse")
  
  ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######
  
  # Comprobar cuántas variantes quedan tras el filtrado
  n_var_chr <- length(gr.filt_chr)
  
  # Guardar nº de variantes por chr
  save(n_var_chr, file = paste0(result_dir, "/n_var_", chr_name, ".RData"))
  # Información en consola:
  cat(paste0("Número de variantes en ", chr_name, ": ", length(gr.filt_chr), "\n"))
  
  ############################################
  ### 4.2. Filtrar los BAMs para sacar la información del cromosoma estudiado
  
  # Establecer el área de búsqueda para el cromosoma concreto
  searchArea_chr <- GRanges(seqnames = c(chr),
                              ranges = IRanges(start = 1, end = (chr_length)))
  
  # Cargar la información del RNAseq del cromosoma estudiado (reads = lecturas de los BAMs)
  reads_chr <- impBamGAL(pathToFiles, searchArea_chr, verbose = TRUE)
  cat("Lecturas del cromosoma ", chr_name, " cargadas exitosamente.\n")
  
  # Mantener solo los primeros 22 cromosomas en reads_chr
  reads_chr <- keepSeqlevels(reads_chr, seqlevels(reads_chr)[1:22], pruning.mode = "coarse")
  # ajustamos la longitud de los cromosomas para que coincidan entre bams y vcf
  seqlengths(gr.filt_chr)[1:22] <- seqlengths(reads_chr)
  
  ############################################
  ### 4.3. Recuento de los alelos que hay en las reads para las posiciones filtradas 
  countList_chr <- getAlleleCounts(reads_chr, gr.filt_chr, verbose = TRUE)
  cat("Recuento de alelos realizado para el cromosoma ", chr_name, ".\n")
  
  # Mostrar el resultado (solo una parte inicial con head)
  print(head(countList_chr))
  
  # Guardar los recuentos resultantes (iniciales o crudos)
  save(countList_chr, file = paste0(result_dir, "/countList_", chr_name, ".RData"))
  cat("Recuendos alélicos del cromosoma ", chr_name, " guardados con éxito .\n")
  
  # Eliminamos las lecturas de los bams para reducir la carga en memoria
  rm(reads_chr)
  cat("Lecturas de los bams eliminadas del enviroment con éxito.")
  
  ############################################
  ### 4.4. Filtrado de los recuentos a partir de la countList
  
  # Inicializa una lista para almacenar los resultados filtrados
  filtered_countList_chr <- list()
  
  # Itera sobre cada posición en countList
  for (pos in names(countList_chr)) {
    # Extrae los conteos de la posición actual
    counts <- countList_chr[[pos]]
    
    # Itera sobre cada fila (cada persona) para verificar los totales
    for (i in 1:nrow(counts)) {
      total_counts <- sum(counts[i, ])  # Suma de lecturas para esa persona
      
      # Si el total es menor que 10, establece toda la fila a 0
      if (total_counts < 10) {
        counts[i, ] <- 0
        cat("Fila ", i, " establecida a 0 para variante ", pos, "\n")
      }
    }
    
    filtered_countList_chr[[pos]] <- counts
    # }
  }
  
  # Mostrar el resultado del filtrado de los recuentos
  cat("Filtrado de recuentos completado para el cromosoma ", chr, ".\n")
  print(head(filtered_countList_chr))
  
  ############################################
  ### 4.5. Crear el objeto ASEset inicial a partir de la lista de conteos filtrada
  
  ase_chr <- ASEsetFromCountList(rowRanges = gr.filt_chr,
                                   filtered_countList_chr)
  cat("Objeto ASEset inicial creado para el cromosoma ", chr,".\n")
  
  ############################################
  ### 4.6. Añadir información de fase (genotipo: homocigotos o heterocigotos)
        # Info de GT se ha extraido previamente para todas las posiciones del VCF con:
        # geno_fullVCF <- geno(VCF)$GT
  
  # Extraer las posiciones de gr.filt_chr
  positions_chr <- start(gr.filt_chr)
  
  # Filtrar el VCF para las variantes del cromosoma usando las posiciones extraídas
  filt_indices_chr <- which(seqnames(VCF) == chr &
                                start(granges(VCF)) %in% positions_chr)
  
  # Extraer la información de genotipos solo para las variantes filtradas
  geno_chr_filt <- geno_fullVCF[filt_indices_chr,]
  
  # Ajustar la información de fase a las dimensiones de `ase_chr` automáticamente
  required_rows <- nrow(ase_chr)  # Número de filas necesario basado en `ase_chr`
  geno_chr_filt_adjusted <- geno_chr_filt[1:required_rows, ]  # Ajustar filas
  
  # Confirmar dimensiones antes de asignar
  cat("Dimensiones de geno_chr_filt_adjusted:", dim(geno_chr_filt_adjusted), "\n")
  cat("Dimensiones de ase_chr:", dim(ase_chr), "\n")
  
  # Asignar la información de fase al objeto ASEset inicial
  phase(ase_chr) <- geno_chr_filt_adjusted
  
  # Checkear el funcionamiento con:
  cat("Información de fase añadida a ASEset inicial para el cromosoma ", chr, ".\n")
  print(head(phase(ase_chr)))
  
  ############################################
  ### 4.7. Filtrado de Heterocigotos y Recuentos = 0
          # Generación del objeto ASEset definitivo
  
  # Fase (phase information) y recuentos (allele counts)
  phase_info_chr <- phase(ase_chr)
  
  # 1. Primer bucle: Poner a 0 los recuentos si son homocigotos o NA
  # Iterar sobre cada variante en la lista de recuentos
  for (variante in names(filtered_countList_chr)) {
    
    # Recorrer cada individuo en la fase y recuentos
    for (individuo in colnames(phase_info_chr)) {
      
      # Si el individuo es homocigoto ("0/0", "1/1", "0|0", "1|1") o NA
      if (phase_info_chr[variante, individuo] %in% c("0/0", "1/1", "0|0", "1|1", "./.")
          || is.na(phase_info_chr[variante, individuo])) {
        # Poner a 0 todos los recuentos para esa variante y ese individuo
        filtered_countList_chr[[variante]][individuo, ] <- 0
      }
      
      # En caso contrario, es decir, si el individuo es heterocigoto
      # ("0/1", "1/0", "0|1", "1|0"), no hacemos nada (los recuentos se mantienen)
    }
  cat(chr_name, "Variante:", variante, "Recuentos de homocigotos = 0.\n")
  }
  
  # 2. Segundo bucle: Verificar si todos los recuentos son 0 y eliminar variantes si es necesario
  
  # Crear una lista para almacenar las variantes eliminadas y conservadas
  var_delete_chr <- c()
  filtered_countList_chr_def <- list()
  
  for (variante in names(filtered_countList_chr)) {
    
    # Extraer los recuentos para la variante actual
    counts <- filtered_countList_chr[[variante]]
    
    # Verificamos si todos los valores son 0 para todas las personas
    if (all(counts == 0)) {
      # Si todos son ceros, añadimos esta variante a la lista de eliminadas
      var_delete_chr <- c(var_delete_chr, variante)
    } else {
      # Si no todos son ceros, añade esta posición a la lista filtrada
      filtered_countList_chr_def[[variante]] <- counts
    }
  cat("Eliminada la variante ",variante,"por recuentos = 0 del", chr_name,".\n")
  }
  
  # Extraer solo las posiciones numéricas de las variantes eliminadas
  # Asumiendo que los nombres son del formato "chrn_posicion"
  pos_delete_chr <- as.numeric(sub(".*_", "", var_delete_chr))
  
  
  # Filtrar el objeto GRanges eliminando las variantes en las posiciones indicadas
  gr.filt_chr_def <- gr.filt_chr[!start(gr.filt_chr)
                                     %in% pos_delete_chr]
  
  
  # Guardar los recuentos resutlantes
  save(filtered_countList_chr, file = paste0(result_dir, "/filtered_countList_chr"
                                             , chr, ".RData"))
  cat("Recuentos filtrados del cromosoma", chr, " guardados correctamente.\n")
  
  ############################################
  ### 4.8. Crear el objeto ASEset definitivo
  # (con los recuentos filtrados previamente, en dos pasos anteriores)
  
  ase_chr_def <- ASEsetFromCountList(rowRanges = gr.filt_chr_def,
                                   filtered_countList_chr_def)
  cat("Objeto ASEset creado para el ", chr_name, ".\n")
  
  # Añadir la información de fase al objeto ASEset definitivo
  # Solo miramos las variantes que han pasado los filtros previos
  positions_def_chr <- start(gr.filt_chr_def)
  
  # Filtrar el VCF para las variantes del cromosoma usando las posiciones extraídas
  phase_indices_chr <- which(seqnames(VCF) == chr &
                                start(granges(VCF)) %in% positions_def_chr)
  
  # Extraer la información de genotipos solo para las variantes filtradas
  phase(ase_chr_def) <- geno_fullVCF[phase_indices_chr, ]
  
  ############################################
  ### 4.9. INFORMACIÓN FENOTÍPICA:
  
  # Cargar el archivo Excel en R y transformarlo a Data Frame
  feno_data <- DataFrame(Condition = fenotipo$Condition)
  
  # Añadir el identificador de las muestras al Data Frame anterior
  rownames(feno_data) <- fenotipo$Sample
  
  # Incluir el fenotipo en el objeto ASEset
  colData(ase_chr_def) <- feno_data
  
  # Para poder checkear que se ha cargado bien podemos utilizar:
  cat("Información fenotípica añadida con éxito para el", chr_name,"\n")
  print(head(colData(ase_chr_def)))
  
  # O también para un individuo concreto:
  # colData(ase_chr_def)["identificador_individuo.bam", "Condition"]
  
  ############################################
  ### 4.10. INFO REF Y ALT ALLELES
  
  # Extraer los alelos de referencia desde la columna REF del gr.filt correspondiente
  ref_alleles_chr <- as.character(mcols(gr.filt_chr_def)[, "REF"])
  
  # Extraer los alelos alternativos desde la columna ALT del gr.filt correspondiente
  alt_alleles_chr <- sapply(mcols(gr.filt_chr_def)$ALT, function(x) as.character(x))
  
  # Se extraen los alelos REF y ALT de forma distinta porque se almacenan de forma distinta en gr.filt:
  # DNAStringSet (REF): Se puede convertir directamente a caracteres.
  # DNAStringSetList (ALT): Se debe extraer cada elemento antes de convertir.
  
  # Verificar los alelos de referencia extraídos
  print(head(ref_alleles_chr))
  print(head(alt_alleles_chr))
  
  # Asignar estos alelos de referencia y alternativo objeto ASEset
  ref(ase_chr_def) <- ref_alleles_chr
  alt(ase_chr_def) <- alt_alleles_chr
  cat("Información Ref y Alt alleles añadida con éxito para el", chr_name, ".\n")
  
  # Verificar la asignación
  cat("Alelo de referencia y alternativo asignado en ASEset:\n")
  print(head(ref(ase_chr_def)))
  print(head(alt(ase_chr_def)))
  
  ############################################
  ### 4.11. INFO GENOTYPE
  
  # Inferir los genotipos y añadirlos en el campo correspondiente del objeto ASEset
  genotype(ase_chr_def) <- inferGenotypes(ase_chr_def)
  cat("Información genotípica añadida con éxito para el", chr_name, ".\n")
  
  ############################################
  ### 4.12. Check Mapping Bias (average reference allele fraction)
  
  # Sacar la fracción del alelo de referencia
  refFrac_chr <- fraction(ase_chr_def, top.fraction.criteria="ref")
  
  # obtener la media
  refFrac_mean_chr <- mean(refFrac_chr, na.rm = TRUE)
  
  # Imprimir el resultado:
  cat("Fracción media del alelo de referencia:", refFrac_mean_chr, "\n")
  save(refFrac_mean_chr, file = paste0(result_dir, "/refFrac_mean_", chr_name, ".RData"))
  
  ############################################
  ### 4.13. Análisis estadístico y Filtrado por significación estadística
  
  # guardar el objeto ASE set definitivo antes de continuar
  save(ase_chr_def, file = paste0(result_dir, "/ase_chr_def_", chr, ".RData"))
  
  ######################     ######################     ######################   
      ##### Test Binomial (AllelicImbalance-Bioconductor):
  ######################     ######################     ###################### 
  
  binom_results_chr <- binom.test(ase_chr_def, n = '*')
  cat("Test binomial realizado con éxito para el ", chr_name, ".\n")
  
  # Creamos un objeto para almacenar los resutlados significativos
  # Cambiamos NaN por NA para facilitar el manejo de esta información
  binom_results_chr[binom_results_chr == "NaN"] <- NA
  
  # Convierte el resultado a un dataframe
  binom_df_chr <- as.data.frame(binom_results_chr)
  
  # Crea un dataframe vacío para almacenar los resultados significativos
  signif_df_chr <- data.frame(matrix(ncol = ncol(binom_df_chr), nrow = 0))
  
  colnames(signif_df_chr) <- colnames(binom_df_chr)
  
  # Inicializar una lista para almacenar las variantes significativas
  binom_sig_variants_chr <- c()
  
  # Bucle para analizar las filas y columnas
  for (i in 1:nrow(binom_df_chr)) {
    for (j in 1:ncol(binom_df_chr)) {
      if (!is.na(binom_df_chr[i, j]) && binom_df_chr[i, j] <= 0.05) {
        # Almacena el valor significativo en el nuevo dataframe
        signif_df_chr <- rbind(signif_df_chr, data.frame(Sample = rownames(binom_df_chr)[i],
                                                           Variant = colnames(binom_df_chr)[j],
                                                           Value = binom_df_chr[i, j],
                                                           Test = "Binomial"))
        # Añadir el nombre de la variante significativa a la lista
        binom_sig_variants_chr <- c(binom_sig_variants_chr, colnames(binom_df_chr)[j])
      }
    }
  }
  cat("Variantes significativas extraídas con éxito para el", chr_name, ".\n")
  
  # Verifica si no se han encontrado resultados significativos
  if (nrow(signif_df_chr) == 0) {
    # Si no hay resultados significativos, dejar el dataframe vacío
    signif_df_chr <- signif_df_chr(Sample = character(),
                                     Variant = character(),
                                     Value = numeric(),
                                     Test = character(),
                                     Strand = character(),
                                     stringsAsFactors = FALSE)
    
    cat("No hay variantes significativas en el test binomial del ", chr_name, ".\n")
  }
  
  # Imprimir el dataframe de resultados significativos
  print(head(signif_df_chr))
  
  # Guardar los resultados del análsis estadístico
  save(signif_df_chr, file = paste0(result_dir, "/signif_df_", chr_name, ".RData"))
  
  # Guardar el número de variantes significativas con este test estadístico
  num_var_binom_chr <- nrow(signif_df_chr)
  save(num_var_binom_chr, file = paste0(result_dir, "/num_var_binom_", chr_name, ".RData"))
  
  # Guardar los nombres de las variantes únicas significativas del test estadístico
  unique_binom_sig_var_chr <- unique(binom_sig_variants_chr)
  save(unique_binom_sig_var_chr, file = paste0(result_dir, "/unique_binom_sig_var_", chr_name, ".RData"))
  
  cat("Información del test binomial del ", chr_name," guardado con éxito.\n")
  
  #########################################################################################
  #########################################################################################
  # 4.14. IDENTIFICADOR DE LAS VARIANTES CONOCIDAS (identificador rs conocido)
  
  ### SOLUCIÓN: solo con binom.test:
  
  # Filtrar el objeto rowRanges del ASEset para quedarnos solo con las variantes comunes significativas
  filtered_gr_chr <- rowRanges(ase_chr_def)[names(rowRanges(ase_chr_def)) 
                                                %in% unique_binom_sig_var_chr]
  
  # Confirmamos que las longitudes de los cromosomas coincidan:
  seqlengths(filtered_gr_chr)[1:22] <- seqlengths(SNPlocs.Hsapiens.dbSNP144.GRCh38)[1:22]
  
  # Cambiamos las localizaciones por el rs 
  rs_signif_var_chr <- getSnpIdFromLocation(filtered_gr_chr, SNPlocs.Hsapiens.dbSNP144.GRCh38)
  
  # Creo un nuevo objeto ASEset solo con la información de las posiciones significativas
  
    # Identificar las posiciones que están tanto en ase_chr_def como en filtered_gr_chr
  pos_to_keep_chr <- which((rownames(ase_chr_def)) %in% unique_binom_sig_var_chr)
  
  # Filtrar el ASEset original para quedarse solo con las posiciones de interés
  ase_chr_def_v2 <- ase_chr_def[pos_to_keep_chr, ]
  
  # Podemos añadir estos rs al objeto ASEset para utilziarlos en los gráficos posteriores
  rowRanges(ase_chr_def_v2) <- rs_signif_var_chr
  
  # Extraer los rs conocidos en una lista
  rs_known_list_chr <- names(rs_signif_var_chr[!grepl(chr_name,"_", names(rs_signif_var_chr))])
  
  # Guardamos el objeto resultante
  save(rs_known_list_chr, file = paste0(result_dir, "/rs_known_list_", chr_name,".RData"))
  cat("Identificadores rs del cromosoma", chr, "guardados con éxito.\n")
  
  #########################################################################################
  #########################################################################################

  cat(" ##########################################\n",
      "Análisis completado para ", chr_name,"\n",
      "##########################################\n")
  
  # Limpiar el entorno dejando solo los archivos necesarios para la siguiente iteración
  rm(list = setdiff(ls(), c("gr.filt", "geno_fullVCF", "pathToFiles",
                            "fenotipo", "chr_sizes", "VCF", "chr")))
  cat("Limpieza del enviroment del cromosoma", chr, " realizada con éxito.\n")
  
}

cat("Éxito: el análisis ha terminado correctamente.")
#########################################################################################
#########################################################################################
#########################################################################################