# Script para procesar los resultados obtenidos y sacar tablas informativas

# Cargar las librerías necesarias
library(dplyr)
library(biomaRt)

################################################################################
################################################################################

# TABLA 1: Número de variantes por cromosoma

################################################################################
### 1. Crear la tabla vacía donde se van a almacenar los resultados:

summary_table <- data.frame(
  Cromosoma = c("1_1", "1_2", "2_1", "2_2", 3:22),
  Variantes_totales = integer(24),
  Variantes_tras_filtrar = integer(24),
  Variantes_signif_binom_test = integer(24),
  RefFrac_mean = numeric(24)
)

################################################################################
### 2. Extraer la información relevante de los objetos correspondientes:

# Establecer el directorio principal donde están las carpetas de los cromosomas
main_dir <- "C:/Users/inesu/Desktop/TFM/RESULTADOS"

for (chr in c("1_1", "1_2", "2_1", "2_2", 3:22)) {
  
  # Construir el nombre de la subcarpeta para el cromosoma actual
  subfolder <- file.path(main_dir, paste0("resultados_chr", chr))
  
  # Cambiar al directorio de la subcarpeta
  setwd(subfolder)
  
  # Determinar el número del cromosoma base (sin el sufijo de mitad)
  chr_num <- ifelse(chr %in% c("1_1", "1_2"), "1", 
                    ifelse(chr %in% c("2_1", "2_2"), "2", as.character(chr)))
  
  # Generar los nombres de los archivos a cargar
  countList_file <- paste0("countList_chr", chr_num, ".RData")
  ase_chr_def_file <- paste0("ase_chr_def_", chr_num, ".RData")
  unique_binom_sig_var_file <- paste0("unique_binom_sig_var_chr", chr_num, ".RData")
  refFrac_mean_file <- paste0("refFrac_mean_chr", chr_num, ".RData")
  
  # Cargar los archivos necesarios
  load(countList_file)            # Carga countList
  load(ase_chr_def_file)          # Carga ase_chr_def
  load(unique_binom_sig_var_file) # Carga unique_binom_sig_var_chr
  load(refFrac_mean_file)         # Carga refFrac_mean_chr
  
  # Asignar valores a la tabla para el cromosoma actual
  summary_table$Variantes_totales[summary_table$Cromosoma == chr] <- length(countList_chr)
  summary_table$Variantes_tras_filtrar[summary_table$Cromosoma == chr] <- length(seqnames(rowRanges(ase_chr_def)))
  summary_table$Variantes_signif_binom_test[summary_table$Cromosoma == chr] <- length(unique_binom_sig_var_chr)
  summary_table$RefFrac_mean[summary_table$Cromosoma == chr] <- refFrac_mean_chr
  
  # Volver al directorio principal
  setwd(main_dir)
}
# Imprimir la tabla completa por mitades
print(summary_table)

################################################################################
### 3. Generar una tabla con los resultados por cromosomas del 1 al 22:
# Crear una copia de la tabla y combinar las mitades para los cromosomas 1 y 2
summary_table_combined <- summary_table %>%
  mutate(Cromosoma = case_when(
    Cromosoma %in% c("1_1", "1_2") ~ "1",  # Combinar las filas de 1_1 y 1_2 en "1"
    Cromosoma %in% c("2_1", "2_2") ~ "2",  # Combinar las filas de 2_1 y 2_2 en "2"
    TRUE ~ as.character(Cromosoma)
  ))

# Agrupar por Cromosoma y calcular la suma y la media de las fracciones de alelo de referencia
summary_table_combined <- summary_table_combined %>%
  group_by(Cromosoma) %>%
  reframe(
    Variantes_totales = sum(Variantes_totales),
    Variantes_tras_filtrar = sum(Variantes_tras_filtrar),
    Variantes_signif_binom_test = sum(Variantes_signif_binom_test),
    RefFrac_mean = mean(RefFrac_mean),
  ) %>%
  arrange(as.numeric(Cromosoma))

# Imprimir la tabla combinada
print(summary_table_combined)


################################################################################
### 4. Guardar las tablas:

# Guardar la tabla completa como archivo CSV
write.csv(summary_table, "summary_table.csv", row.names = FALSE)

# Guardar la tabla combinada como archivo CSV
write.csv(summary_table_combined, "summary_table_def.csv", row.names = FALSE)

################################################################################
################################################################################
 
# TABLA 2: Genes significativos en el Test de Wilcoxon

################################################################################

# Paso 1: Configurar la conexión con biomaRt
ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", 
                   host = "https://www.ensembl.org")

# Paso 2: Cargar los genes que han salido significativos (excel con gen y cromosoma)
genes <- read_excel("C:/Users/inesu/Desktop/TFM/RESULTADOS/wilcox_genes.xlsx")


# Crear una tabla vacía para la salida final
final_table <- data.frame(
  ensembl_gene_id = character(), 
  external_gene_name = character(),
  chromosoma = character(), 
  start_position = integer(),
  end_position = integer(),
  wilcoxon_pval = numeric(),
  stringsAsFactors = FALSE
)

# Iterar sobre los genes
for (i in 3:nrow(genes)) {
  gene_id <- genes$ID[i]
  chr <- genes$CHR[i]
  
  # Cargar los resultados de Wilcoxon para el cromosoma correspondiente
  subfolder <- paste0("C:/Users/inesu/Desktop/TFM/RESULTADOS/resultados_chr", chr)
  wilcox_result_file <- paste0(subfolder, "/wilcox_results_chr", gsub("_", "", chr), ".RData")
  
  # Cargar el archivo de resultados de Wilcoxon
  load(wilcox_result_file)
  
  # Buscar el p-valor del gen en los resultados de Wilcoxon
  p_value <- wilcox_results$p_value[wilcox_results$ensembl_gene_id == gene_id]
  
  # Obtener la información del gen desde Ensembl
  gene_info <- getBM(
    attributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name", "start_position", "end_position"),
    filters = "ensembl_gene_id",
    values = gene_id,
    mart = ensembl
  )
  
  # Si no se encuentra la información del gen, continuar con el siguiente gen
  if (nrow(gene_info) == 0) {
    next
  }
  
  # Crear una fila para agregar a la tabla final
  new_row <- data.frame(
    ensembl_gene_id = gene_info$ensembl_gene_id[1],
    external_gene_name = gene_info$external_gene_name[1],
    chromosoma = gene_info$chromosome_name[1],
    start_position = gene_info$start_position[1],
    end_position = gene_info$end_position[1],
    wilcoxon_pval = p_value,
    stringsAsFactors = FALSE
  )
  
  # Añadir la fila a la tabla final
  final_table <- rbind(final_table, new_row)
}

# Ver las primeras filas de la tabla final
print(head(final_table))

# Guardar la tabla final como CSV
write.csv(final_table, "tabla_final_genes_wilcoxon.csv", row.names = FALSE)
