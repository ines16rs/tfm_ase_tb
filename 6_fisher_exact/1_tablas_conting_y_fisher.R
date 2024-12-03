
################################################################################
# Cargar AllelicImbalance, VariantAnnotation y otros paquetes necesarios
library(AllelicImbalance)
library(SNPlocs.Hsapiens.dbSNP144.GRCh38)
library(biomaRt)
library(dplyr)

################################################################################
# Establecer la ruta donde queremos guardar los resultados
result_dir <- getwd()
data_dir <- "C:/Users/inesu/Desktop/TFM/wilcoxon/data/"
subcarpeta <- "C:/Users/inesu/Desktop/TFM/wilcoxon/recuentos_concretos/"

# Crear una tabla de fenotipos a partir de colData
fenotipos <- as.data.frame(colData(ase_chr_def))

# Bucle que recorre cada valor en el vector (los cromosomas estudiados)
for (valor in c(22, 21, 19:3, "2_1", "2_2", "1_2", "1_1")) {
  
  # Notificación del inicio del proceso para cada cromosoma
  cat("Trabajando con el cromosoma:", valor, "\n")
  
  # Cargar el archivo .RData correspondiente
  load(file.path(data_dir, paste0("ase_chr_def_", valor, ".RData")))
  cat("Archivo ase_chr_def_",valor, "cargado correctamente.\n")
  
  setwd(data_dir)
  # Cargar el archivo .txt correspondiente con las variantes de interés
  # Generar el nombre del archivo dinámicamente
  file_name <- paste0("signif_var_wilcox_chr", valor, ".txt")
  signif_data <- read.table(file_name, header = TRUE, sep = "\t")
  var_list <- list()
  var_153 <- var_list[[paste0("chr", valor)]] <- paste0("chr", signif_data$Cromosoma, "_", signif_data$Posicion)
  cat("Archivo .txt del cromosoma", valor, " cargado correctamente.\n")
  setwd(result_dir)
  
  ################################################################################  
  # Crear un data frame vacío donde se guardarán los resultados
  proporciones_medias <- data.frame(Cromosoma = character(),
                                    Posicion = integer(),
                                    Affected = numeric(),
                                    Unaffected = numeric(),
                                    cambio_de_alelo = character(),
                                    stringsAsFactors = FALSE)
  ################################################################################  
  
  # Iterar sobre las variantes significativas
  for (variant in var_153) {
    
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
    
    # Obtener los recuentos alélicos
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
    
    ############################################################################
    
    # Crear una nueva columna para clasificar el alelo mayoritario
    datos_variante$AleloMayoritario <- ifelse(datos_variante$ProporcionRef > 0.5, 
                                              "Referencia", "Alternativo")
    print(datos_variante)
    
    # Crear la tabla de contingencia
    tabla_contingencia <- table(datos_variante$Condition, datos_variante$AleloMayoritario)
    print(tabla_contingencia)
    
    # Verificar las dimensiones de la tabla
    dim_tabla <- dim(tabla_contingencia)
    
    # Si la tabla no tiene 2 filas, añadir una fila de ceros
    if (dim_tabla[1] < 2) {
      tabla_contingencia <- rbind(tabla_contingencia, c(0, 0))  # Añadir una fila de ceros
    }
    
    # Si la tabla no tiene 2 columnas, añadir una columna de ceros
    if (dim_tabla[2] < 2) {
      tabla_contingencia <- cbind(tabla_contingencia, c(0, 0))  # Añadir una columna de ceros
    }
    
    write.table(
      tabla_contingencia,
      file = paste0("tabla_conting_", variant, ".txt"),
      sep = "\t",
      row.names = TRUE,
      col.names = NA,
      quote = FALSE
    )
    
    # Realizar el test exacto de fisher (es un tipo de test de asociación)
    # Utilizamos Fisher porque es como Chi cuadrado, pero para nº datos reducido
    fisher_result <- fisher.test(tabla_contingencia)
    print(fisher_result)
    
    # Extraer p-valor y odds ratio del test de Fisher
    p_valor_fisher <- fisher_result$p.value
    odds_ratio <- ifelse(is.null(fisher_result$estimate), NA, fisher_result$estimate)
    ############################################################################
    
    medias_por_grupo <- datos_variante %>%
      group_by(Condition) %>%
      summarize(media = mean(ProporcionRef, na.rm = TRUE))
    
    media_df = as.data.frame(t(medias_por_grupo$media))
    
    # Asignar nombres a las columnas para facilitar su comprensión
    colnames(media_df) <- c("Affected", "Unaffected")
    
    # Añadir la columna cromosoma, posición, p-val y odd ratio de la variante al resultado
    media_df$Cromosoma <- cromosoma_var
    media_df$Posicion <- posicion_var
    media_df$Pvalor_Fisher <- p_valor_fisher
    media_df$Odds_Ratio <- odds_ratio
    
    # Verificar el cambio de la expresión predominante del alelo (SI O NO)
    media_df$cambio_de_alelo <- ifelse(media_df$Pvalor_Fisher <= 0.05, "SI", "NO")
    
    
    # Guardar allele_counts_muestras en un archivo .txt solo si cambio_de_alelo es "SI"
    if (media_df$cambio_de_alelo == "SI") {
      
      setwd(subcarpeta)
      
      # Guardar el objeto allele_counts_muestras en un archivo .txt
      write.table(
        allele_counts_muestras,
        file = paste0("subset_counts_", variant, ".txt"),
        sep = "\t",
        row.names = TRUE,
        col.names = NA,
        quote = FALSE
      )
      setwd(result_dir)
    }
    
    # Reorganizar las columnas para que la tabla tenga la estructura correcta
    media_df <- media_df[, c("Cromosoma", "Posicion", "Affected", "Unaffected",
                             "Pvalor_Fisher", "cambio_de_alelo", "Odds_Ratio")]
    
    # Añadir el resultado de la variante a la tabla final
    proporciones_medias <- bind_rows(proporciones_medias, media_df)
    
  }  
  
  # Guardar las tablas solo si contienen datos
  if (nrow(proporciones_medias) > 0) {
    # Crear el nombre del archivo con la iteración actual (valor)
    file_name <- paste0("data_solo_fisher_", valor, ".txt")
    write.table(proporciones_medias, file = file_name, sep = "\t", row.names = FALSE, quote = FALSE)
    cat("Se ha generado la tabla de medias para el cromosoma:", valor, "\n")
  }
  
  head(proporciones_medias)
  
  # Realizar la unión de signif_data con proporciones_medias
  signif_data$Cromosoma <- as.character(signif_data$Cromosoma)
  proporciones_medias$Cromosoma <- as.character(proporciones_medias$Cromosoma)
  # Añadir columnas solo cuando las columnas "Cromosoma" y "Posicion" coinciden
  signif_data <- signif_data %>%
    left_join(proporciones_medias %>% select(Cromosoma, Posicion, Affected,
                                             Unaffected, Pvalor_Fisher,
                                             cambio_de_alelo, Odds_Ratio),
              by = c("Cromosoma", "Posicion"))
  
  tabla_final <- paste0("fisher_full_", valor, ".txt")
  write.table(signif_data, file = tabla_final, sep = "\t", row.names = FALSE, quote = FALSE)
  cat("Se ha generado la tabla de final para el cromosoma:", valor, "\n")
  
  # Eliminamos todos los objetos que no sean necesarios
  rm(list = setdiff(ls(), c("fenotipos", "result_dir", "data_dir",
                            "subcarpeta", "ase_chr_def", "valor", "proporciones_medias")))
  
}