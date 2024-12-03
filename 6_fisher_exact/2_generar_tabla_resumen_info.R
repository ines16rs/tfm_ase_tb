
library(dplyr)
library(writexl)

# Definir el directorio donde están los archivos .txt de partida
dir_path <- "C:/Users/inesu/Desktop/TFM/wilcoxon/Resultados/full_data"

# Obtener una lista de todos los archivos que coinciden con el patrón
files <- list.files(path = dir_path, pattern = "^medias_wilcox_.*\\.txt$", full.names = TRUE)

# Leer el primer archivo con encabezado
first_data <- read.table(files[1], header = TRUE, sep = "\t")

# Leer los archivos restantes sin encabezado, asignando los nombres de columna del primer archivo
other_data <- lapply(files[-1], function(file) {
  data <- read.table(file, header = FALSE, sep = "\t")
  colnames(data) <- colnames(first_data)
  return(data)
})

# Combinar el primer archivo con los demás archivos en un solo data frame
will_wilcoxon <- do.call(rbind, c(list(first_data), other_data))

# Eliminar las filas adicionales del encabezado, si existen
will_wilcoxon <- will_wilcoxon[will_wilcoxon$rsID != "rsID", ]


# Guardar el resultado en un archivo Excel y en otro RData
write_xlsx(will_wilcoxon, "will_wilcoxon.xlsx")
save(will_wilcoxon, file = "will_wilcoxon.RData")
