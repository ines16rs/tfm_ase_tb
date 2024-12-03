setwd("C:/irs_tfm")
getwd()

################################################################################
# Cargar AllelicImbalance, VariantAnnotation y otros paquetes necesarios
library(AllelicImbalance)
library(SNPlocs.Hsapiens.dbSNP144.GRCh38)
library(biomaRt)
library(dplyr)

################################################################################
# Establecer la ruta donde vamos a guardar los resultados de este paso del TFM
result_dir <- getwd()
data_dir <- "C:/Users/inesu/Desktop/TFM/wilcoxon/data/"

# Conectar con Ensembl para obtener información de genes
ensembl_gen <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", 
                   host = "https://www.ensembl.org")

# Crear una tabla de fenotipos a partir de colData
fenotipos <- as.data.frame(colData(ase_chr_def))

# Crear una tabla para almacenar el conteo de variantes en cada cromosoma
tabla_conteo <- data.frame(Cromosoma = character(),
                           N_var_analizadas = integer(),
                           N_var_signif_wilcoxon = integer(),
                           stringsAsFactors = FALSE)

# Bucle que recorre cada valor en el vector (son los cromosomas)
# En el caso del cromosoma 1 y 2, como eran muy grandes se han hecho por mitades
for (valor in c("2_1", "2_2", "1_2", "1_1", 3:22)) {
  
  # Notificación del inicio del proceso para cada cromosoma
  cat("Trabajando con el cromosoma:", valor, "\n")
  
  # Cargar el archivo .RData correspondiente
  load(file.path(data_dir, paste0("ase_chr_def_", valor, ".RData")))
  cat("Archivo ase_chr_def_", valor, " cargado correctamente.\n")
  
  # Cargar el archivo .RData correspondiente
  load(file.path(data_dir, paste0("signif_df_chr", valor, ".RData")))
  cat("Archivo signif_df_chr", valor, " cargado correctamente.\n")
  
  ################################################################################  
  # Data frame para almacenar los resultados finales
  tabla_resultados <- data.frame(rsID = character(),
                                 Cromosoma = character(),
                                 Posicion = integer(),
                                 pvalor_wilcox = numeric(),
                                 pvalor_fdr = numeric(),
                                 Gen_Ensembl_ID = character(),
                                 Gen_External_Name = character(),
                                 Num_individuos = character(),
                                 stringsAsFactors = FALSE)
  
  # Reestablecer el contador de variantes analizadas en el test de Wilcoxon
  contador_variantes <- 0
  contador_variantes_signif <- 0
  
  # Iterar sobre las variantes significativas
  for (variant in unique(signif_df_chr$Variant)) {
    
    cat("Analizando variante:", variant, "\n")
    
    # Obtener el índice de la variante en ase_chr_def
    variant_idx <- match(variant, rownames(ase_chr_def))
    if (is.na(variant_idx)) {
      warning(paste("La variante", variant, "no se encontró en ase_chr_def."))
      next
    }
    
    # Extraer el cromosoma y la posición directamente del nombre de la variante
    split_variant <- strsplit(variant, "_")[[1]]
    cromosoma_var <- gsub("chr", "", split_variant[1])  # Parte antes del "_"
    posicion_var <- as.integer(split_variant[2])  # Parte después del "_"
    # Crear la región cromosómica en el formato correcto para BioMart
    region <- paste0(cromosoma_var, ":", posicion_var, "-", posicion_var)
    
    # Extraer los alelos de referencia y alternativos de rowData
    ref_allele <- rowData(ase_chr_def)$ref[variant_idx]
    alt_allele <- rowData(ase_chr_def)$alt[variant_idx]
    
    # Verifica que el alelo de referencia sea A, C, G o T
    if (!(ref_allele %in% c("A", "C", "G", "T"))) {
      # Si no es A, C, G o T, pasa a la siguiente iteración o excluye la variante
      next # usa `next` si estás en un bucle o ajusta según tu código para omitir la variante
    }
    
    # Verifica que el alelo alternativo sea A, C, G o T
    if (!(alt_allele %in% c("A", "C", "G", "T"))) {
      # Si no es A, C, G o T, pasa a la siguiente iteración o excluye la variante
      next # usa `next` si estás en un bucle o ajusta según tu código para omitir la variante
    }
    
    # Obtener los recuentos alélicos solo para las muestras relevantes
    allele_counts <- alleleCounts(ase_chr_def)[[variant]]
    
    # Sumar los recuentos por cada persona (fila)
    suma_recuentos <- rowSums(allele_counts)
    
    # Crear un vector con los identificadores de las personas que tienen al menos un recuento distinto de 0
    muestras_para_variante <- rownames(allele_counts)[suma_recuentos > 0]
    
    allele_counts_muestras <- allele_counts[muestras_para_variante, , drop = FALSE]
    
    # Calcular proporciones de alelo de referencia, manejando errores
    cat("Realizando cálculos de proporciones vara la variante:", variant, "\n")
    total_counts <- allele_counts_muestras[, ref_allele] + allele_counts_muestras[, alt_allele]
    proporciones_ref <- ifelse(total_counts == 0, NA, allele_counts_muestras[, ref_allele] / total_counts)
    
    # Si el alelo de referencia tiene recuentos 0, asignar proporción 0
    proporciones_ref[allele_counts_muestras[, ref_allele] == 0] <- 0
    
    # Extraer los fenotipos de las muestras relevantes a partir de la tabla de fenotipos
    fenotipos_muestras <- fenotipos[muestras_para_variante, "Condition", drop = FALSE]
    
    # Crear un dataframe con las proporciones de alelo de referencia y el fenotipo
    datos_variante <- data.frame(Sample = rownames(fenotipos_muestras),
                                 ProporcionRef = as.vector(proporciones_ref),
                                 Condition = fenotipos_muestras$Condition)
    print(datos_variante)
    
    # Contar el número de muestras por condición
    num_sanos <- sum(datos_variante$Condition == "Unaffected")
    num_enfermos <- sum(datos_variante$Condition == "Affected")
    
    # Solo realizar el test si hay al menos 2 muestras sanas y 2 enfermas
    if (num_sanos >= 2 & num_enfermos >= 2) {
      cat("Realizando test de Wilcoxon para la variante:", variant, "\n")
      
      # Incrementar el contador de variantes
      contador_variantes <- contador_variantes + 1
      
      # Realizar el test de Wilcoxon Rank Sum o U de Mann Whitney
      wilcox_test <- wilcox.test(ProporcionRef ~ Condition, data = datos_variante)
      cat("Test de wilcoxon realizado con éxito.", "\n")
      
      # Guardar siempre los resultados del test, no solo los significativos
      tabla_resultados <- rbind(tabla_resultados, data.frame(
        rsID = NA, Cromosoma = cromosoma_var, 
        Posicion = posicion_var, 
        pvalor_wilcox = wilcox_test$p.value,
        pvalor_fdr = NA, Gen_Ensembl_ID = NA, 
        Gen_External_Name = NA, Num_individuos = paste(num_sanos, num_enfermos, sep = " + ")
      ))
      
      # Verificar si el test es significativo (p-valor < 0.05)
      if (!is.na(wilcox_test$p.value) && wilcox_test$p.value < 0.05) {
        cat("WILCOXON TEST SIGNIFICATIVO para la variante:", variant, "\n")
        cat("Buscando información del SNP y el gen para la variante:", variant, "\n")
        
        contador_variantes_signif <- contador_variantes_signif + 1
        
        # Obtener información adicional sobre la variante desde Ensembl
        variant_info <- getBM(attributes = c("chromosome_name", "start_position", 
                                             "ensembl_gene_id", "external_gene_name"),
                              filters = c("chromosome_name", "start", "end"), 
                              values = list(cromosoma_var,posicion_var,posicion_var),
                              mart = ensembl_gen)
        
        # Genero un objeto gr para buscar los rsIDs con SNPlocs.Hsapiens.dbSNP144.GRCh38
        gr <- GRanges(seqnames = cromosoma_var, ranges = IRanges(start = posicion_var, end = posicion_var))
        # Obtener los rsIDs de las variantes usando getSnpIdFromLocation
        gr_with_rs <- getSnpIdFromLocation(gr,SNPlocs.Hsapiens.dbSNP144.GRCh38)
        # Lo pasamos a un dataframe y el rs quedará como el rowname
        df_rs <- as.data.frame(gr_with_rs)
        
        # Si se encuentra información de la variante, añadir a la tabla de resultados
        if (nrow(variant_info) > 0) {
          
          cat("Añadiendo información del SNP y el gen para la variante:", variant, "\n")
          # Obtener los datos que necesitaremos a continuación
          rsID <- rownames(df_rs)
          cromosoma <- cromosoma_var
          posicion <- posicion_var
          gen_ensembl_id <- variant_info$ensembl_gene_id[1]
          gen_name <- variant_info$external_gene_name[1]
          
          # Formato del número de individuos como "sanos + enfermos"
          num_individuos <- paste(num_sanos, num_enfermos, sep = " + ")
          
          # Buscar la fila correspondiente en tabla_resultados usando la columna 'Posicion'
          fila_index <- which(tabla_resultados$Posicion == posicion_var & tabla_resultados$Cromosoma == cromosoma_var)
          
          # Comprobar si se ha encontrado la fila
          if (length(fila_index) == 1) {
            # Actualizar la fila encontrada con la nueva información
            tabla_resultados[fila_index, "rsID"] <- rsID
            tabla_resultados[fila_index, "Gen_Ensembl_ID"] <- gen_ensembl_id
            tabla_resultados[fila_index, "Gen_External_Name"] <- gen_name
            tabla_resultados[fila_index, "Num_individuos"] <- num_individuos
          } else {
            warning("No se ha encontrado la fila correspondiente en la tabla para la posición ", posicion_var)
          }
        }
      }
    } else {
      # Si no hay suficientes muestras, registrar la variante con solo el p-valor
      cat("No hay suficientes muestras para analizar la variante", variant,".\n")
      cat("Número de sanos =", num_sanos,". Número de enfermos=", num_enfermos,".")
    }
  }  
  
  # Guardar el número de variantes analizadas para este cromosoma
  tabla_conteo <- rbind(tabla_conteo,
                        data.frame(
                          Cromosoma = valor,
                          N_var_analizadas = contador_variantes,
                          N_var_signif_wilcoxon = contador_variantes_signif))
  cat("Tabla de conteos actualizada con:", contador_variantes,"y",contador_variantes_signif, "\n")
  
  # Ajustar los p-valores mediante FDR 
  # (por el peq. tamaño muestral no se esperan resultados significativos)
  tabla_resultados$pvalor_fdr <- p.adjust(tabla_resultados$pvalor_wilcox, method = "fdr")
  archivo_restultados <- paste0("resultados_wilcox_chr", valor, ".txt")
  write.table(tabla_resultados, file = archivo_restultados, sep = "\t", row.names = FALSE, quote = FALSE)
  cat("Se han guardado los resultados obtenidos del cromosoma:", valor, "\n")
  
  # Filtrar las variantes significativas
  tabla_significativa <- tabla_resultados[!is.na(tabla_resultados$pvalor_wilcox) & tabla_resultados$pvalor_wilcox < 0.05, ]
  tabla_significativa_fdr <- tabla_resultados[!is.na(tabla_resultados$pvalor_fdr) & tabla_resultados$pvalor_fdr < 0.05, ]
  
  # Guardar las tablas solo si contienen datos
  if (nrow(tabla_significativa) > 0) {
    # Crear el nombre del archivo con la iteración actual (valor)
    archivo_significativa <- paste0("signif_var_wilcox_chr", valor, ".txt")
    write.table(tabla_significativa, file = archivo_significativa, sep = "\t", row.names = FALSE, quote = FALSE)
    cat("Se ha generado la tabla de variantes significativas para el cromosoma:", valor, "\n")
  }
  
  if (nrow(tabla_significativa_fdr) > 0) {
    # Crear el nombre del archivo con la iteración actual (valor)
    archivo_significativa_fdr <- paste0("signif_var_fdr_chr", valor, ".txt")
    write.table(tabla_significativa_fdr, file = archivo_significativa_fdr, sep = "\t", row.names = FALSE, quote = FALSE)
    cat("Se ha generado la tabla de variantes SIGNIF-FDR para el cromosoma:", valor, "\n")
  }
  
  # Mostrar las tablas de resultados significativos
  head(tabla_resultados)
  cat("Tabla resultados con", dim(tabla_resultados)[1],"variantes analizadas y",
      dim(tabla_resultados)[2], "personas incluidas.", "\n")
  print(tabla_significativa)
  print(tabla_significativa_fdr)
  
  # Eliminamos todos los objetos que no sean necesarios
  rm(list = setdiff(ls(), c("fenotipos", "result_dir", "ensembl_gen", "valor", "data_dir", "tabla_conteo")))
}

# Guardar la tabla de conteo en un archivo
write.table(tabla_conteo, file = "tabla_conteo_variantes.txt", sep = "\t", row.names = FALSE, quote = FALSE)
cat("Se ha generado la tabla de conteo de variantes analizadas.\n")
