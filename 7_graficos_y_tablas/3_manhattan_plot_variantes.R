# Cargar las librerías necesarias
library(ggplot2)
library(dplyr)
library(readr)

data_dir = "C:/irs_tfm/wilcoxon/Resultados_2_sanos_vs_enf/"

# Crear una lista de archivos a cargar
archivos <- list.files(path = data_dir, pattern = "resultados_wilcox_chr.*\\.txt", full.names = TRUE)

# Leer y combinar todos los archivos
datos <- archivos %>%
  lapply(read_delim, delim = "\t") %>%  # Leer cada archivo
  bind_rows(.id = "Archivo") %>%  # Combinar en un solo dataframe
  mutate(Cromosoma = factor(Cromosoma, levels = c(1, "1_1", "1_2", 2, "2_1", "2_2", 3:22))) # Orden correcto

# Agrupar los cromosomas 1_1 y 1_2, así como 2_1 y 2_2, en un solo cromosoma
datos <- datos %>%
  mutate(Cromosoma = case_when(
    Cromosoma %in% c("1_1", "1_2") ~ "1",
    Cromosoma %in% c("2_1", "2_2") ~ "2",
    TRUE ~ as.character(Cromosoma)
  ))

# Suponiendo que tus datos están almacenados en un objeto llamado 'datos'
datos_seleccionados <- datos %>%
  select(rsID, Cromosoma, Posicion, pvalor_wilcox)

# Convertir las columnas a tipo numérico
datos_seleccionados_numeric <- datos_seleccionados %>%
  mutate(
    Cromosoma = as.numeric(Cromosoma),
    Posicion = as.numeric(Posicion),
    pvalor_wilcox = as.numeric(pvalor_wilcox)
  )

################################################################################
# 1. Cada cromosoma por separado: 
#####     #####     #####     #####     #####     #####     #####     #####     
# Crear una columna de color para cada cromosoma (pares en skyblue y impares en otro tono de azul)
datos_seleccionados_numeric <- datos_seleccionados_numeric %>%
  mutate(
    Cromosoma = as.factor(Cromosoma),  # Asegurarse de que Cromosoma sea un factor
    color_cromosoma = ifelse(as.numeric(Cromosoma) %% 2 == 0, "skyblue", "dodgerblue")
  )

# Crear el Manhattan plot
ggplot(datos_seleccionados_numeric, aes(x = Posicion, y = -log10(pvalor_wilcox), color = color_cromosoma)) +
  geom_point(alpha = 0.7) +  # Usar puntos transparentes
  scale_color_identity() +  # Usar colores definidos por el usuario
  facet_wrap(~ Cromosoma, scales = "free_x", ncol = 6) +  # Un gráfico por cromosoma
  theme_minimal() +
  theme(
    strip.text = element_text(size = 8),  # Tamaño de texto para los cromosomas
    axis.text.x = element_text(angle = 90, hjust = 1),  # Rotar etiquetas del eje X
    axis.title.x = element_blank(),  # Quitar título eje X
    axis.title.y = element_text(size = 10),
    legend.position = "none"  # Eliminar leyenda
  ) +
  labs(y = "-log10(p-valor)", title = "Manhattan Plot") +
  theme(
    plot.title = element_text(hjust = 0.5)  # Centrar título
  )

################################################################################
# 2: Todos los cromosomas juntos:
#####     #####     #####     #####     #####     #####     #####     #####     
# Datos de tamaño de los cromosomas
chr_sizes <- c(248956422, 242193529, 198295559, 190214555, 181538259, 
               170805979, 159345973, 145138636, 138394717, 133797422, 
               135086622, 133275309, 114364328, 107043718, 101991189, 
               90338345, 83257441, 80373285, 58617616, 64444167, 
               46709983, 50818468)

# Filtrar datos eliminando filas con NA y valores de pvalor_wilcox fuera del rango [0,1]
datos_seleccionados_numeric <- datos_seleccionados_numeric %>%
  filter(!is.na(pvalor_wilcox) & pvalor_wilcox >= 0 & pvalor_wilcox <= 1) %>%  # Filtrar NA y p-valor fuera de rango
  mutate(
    Cromosoma = as.factor(Cromosoma),  # Asegurarse de que Cromosoma sea un factor
    color_cromosoma = ifelse(as.numeric(Cromosoma) %% 2 == 0, "skyblue", "dodgerblue"),
    Cromosoma_numerico = as.numeric(Cromosoma)  # Crear columna numérica para ordenar cromosomas
  ) %>%
  # Ajustar las posiciones según el tamaño de los cromosomas
  mutate(
    pos_cromosoma = Posicion + cumsum(c(0, chr_sizes[-length(chr_sizes)]))[Cromosoma_numerico]
  )

ggplot(datos_seleccionados_numeric, aes(x = pos_cromosoma, y = -log10(pvalor_wilcox), color = color_cromosoma)) +
  geom_point(alpha = 0.7, size = 1) +  # Usar puntos con transparencia
  scale_color_identity() +  # Usar colores definidos por el usuario
  scale_x_continuous(
    breaks = cumsum(chr_sizes) - chr_sizes / 2,  # Poner un corte en el centro de cada cromosoma
    labels = as.character(1:22),  # Etiquetas para los cromosomas, solo el número
    expand = c(0, 0)  # Eliminar espacio extra antes del primer cromosoma
  ) +
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "red", linewidth = 1.0) +  # Línea de significación más gruesa
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 7, angle = 0, hjust = 0.5, vjust = 2),  # Reducir tamaño y ajustar la posición de las etiquetas
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 10),
    legend.position = "none",  # Eliminar leyenda
    axis.ticks.x = element_line(size = 0.5)  # Añadir algo de grosor a las marcas del eje
  ) +
  labs(
    y = "-log10(p-valor)",
    x = "Cromosoma",  # Título del eje X
    title = "Manhattan Plot - Diferencias en Expresión Alelica"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5)  # Centrar título
  )

# Guardar el gráfico en un archivo PNG de alta calidad
ggsave("manhattan_plot2.png", 
       plot = last_plot(), 
       width = 10, 
       height = 3, 
       dpi = 300)  # Alta resolución (300 dpi)

ggsave("manhattan_plot.pdf", plot = last_plot(), width = 10, height = 4)
################################################################################
################################################################################
################################################################################
