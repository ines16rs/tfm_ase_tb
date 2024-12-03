# Script para procesar los bams obtenidos y sacar el número de reads en cada uno

# Cargar las librerías necesarias
library(readxl)
library(ggplot2)
library(tidyr)
library(dplyr)
library(scales)

################################################################################
################################################################################

### VERSIÓN en ESPAÑOL:

# Cargar el archivo Excel con los datos necesarios (número de lecturas en bams)
# (Estos datos se pueden sacar de los archivos flagstats generados en bash y del
# número de lecturas de los archivos generados en el paso 5_merge de WASP)
data <- read_excel("C:/irs_tfm/RESULTADOS/TABLAS/recuentos_bams_tabla_boxplots.xlsx")
head(data)

# Renombrar las columnas para evitar problemas con caracteres especiales
colnames(data) <- c("Sample", "Condition", "Nº_Lecturas_Iniciales", "Nº_Lecturas_Sin_Duplicados", "Nº_Lecturas_WASP")

# Convertir los datos a formato largo para ggplot
data_long <- data %>%
  pivot_longer(cols = starts_with("Nº_Lecturas"),  # Selecciona las columnas de lecturas
               names_to = "Tipo_lectura",           # El nombre de la columna de tipos de lecturas
               values_to = "Lecturas")              # Los valores correspondientes a las lecturas

# Reemplazar los nombres de las categorías en la columna Tipo_lectura por los nuevos nombres
data_long$Tipo_lectura <- recode(data_long$Tipo_lectura,
                                 "Nº_Lecturas_Iniciales" = "Archivos crudos",
                                 "Nº_Lecturas_Sin_Duplicados" = "Sin duplicados",
                                 "Nº_Lecturas_WASP" = "Post-WASP")

# Reordenar las categorías en 'Tipo_lectura' para que 'Lecturas iniciales' esté primero
data_long$Tipo_lectura <- factor(data_long$Tipo_lectura, levels = c("Archivos crudos", "Sin duplicados", "Post-WASP"))

# Reemplazar los valores de 'Condition' para traducirlos
data_long$Condition <- recode(data_long$Condition,
                              "Affected" = "Afectados",
                              "Unaffected" = "No afectados")

# Crear el boxplot con ggplot2 en español
ggplot(data_long, aes(x = Tipo_lectura, y = Lecturas, fill = Condition)) +
  geom_boxplot() +                         # Tipo de gráfico: boxplot
  theme_minimal() +                        # Tema de fondo
  labs(title = "Boxplots de lecturas por tipo y condición",
       x = "Tipo de archivos", 
       y = "Número de lecturas") +          # Etiquetas en español
  scale_y_continuous(labels = label_comma()) + # Eliminar formato científico en el eje Y
  scale_fill_manual(values = c("Afectados" = "coral", "No afectados" = "skyblue")) + # Colores personalizados
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),  # Sin rotar los nombres del eje X
        legend.title = element_text(size = 12),             # Ajustar tamaño de título de leyenda
        legend.text = element_text(size = 10))              # Ajustar tamaño del texto de la leyenda

# Guardar el gráfico generado
ggsave("boxplots_reads_by_type_spa.png", width = 10, height = 6)
