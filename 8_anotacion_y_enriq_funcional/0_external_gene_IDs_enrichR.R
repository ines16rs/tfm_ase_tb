
getwd()

###############################################################################
###############################################################################
### Cargar el archivo que conteine los IDs estables (Stable Gene ID o ENSXXXXXX)

file_genes <- "genes_ENSG.txt"

table_genes <- read.table(file_genes, header = FALSE, stringsAsFactors = FALSE, sep = "\t")

nrow((table_genes)) # Ver el número de genes que tenemos
# View(table_genes) # Comprobar que se ha cargado bien

################################################
### Instalar y cargar BiomaRt

# install.packages("biomaRt")
library(biomaRt)

## Conectar con Ensembl-Biomart usando hg38 como genoma de referencia humano

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://www.ensembl.org")

################################################
### Convertir los IDs estables a símbolos de genes Entrez utilizando Biomart

## Resultados - Atributos de interés:
## Gene id = "ensembl_gene_id"
## Entrez gene symbol = "external_gene_name"

genes <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
               filters = "ensembl_gene_id",
               values = table_genes,
               mart = ensembl)

nrow((genes)) # Comprobar que la tabla sigue teniendo el mismo nº de genes
# View(genes) # Ver que el formato de los resultados es el correcto

## Guardar los símbolos de genes Entrez en un archivo de texto
write.table(unique(genes$external_gene_name), "entrez_gene_symbols.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

# lista de genes
gene_list <- genes$external_gene_name


###############################################################################
###############################################################################

### ENRIQUECIMIENTO FUNCIONAL - EnrichR web
### https://maayanlab.cloud/Enrichr/enrich

###############################################################################