# Cargar paquetes necesarios
library(edgeR)

# Supongamos que tienes un archivo CSV llamado 'datos_miRNA.csv' con los datos de expresión
# Asegúrate de que los datos estén en el formato correcto para edgeR, por ejemplo, TPM o FPKM

# Cargar los datos de expresión de microARNs
datos_expresion <- read.csv("datos_miRNA.csv", header=TRUE, row.names=1)

# Crear el objeto DGEList para edgeR
dge <- DGEList(counts=datos_expresion)

# Filtrar y normalizar los datos
dge <- calcNormFactors(dge)
dge <- estimateCommonDisp(dge)
dge <- estimateTagwiseDisp(dge)

# Diseño experimental (ejemplo: comparación entre dos condiciones)
grupo <- factor(c(rep("Condición1", 3), rep("Condición2", 3)))
design <- model.matrix(~grupo)

# Ajustar el modelo y realizar el análisis diferencial
fit <- glmFit(dge, design)
contrast <- makeContrasts(grupoCondición2 - grupoCondición1, levels=design)
res <- glmLRT(fit, contrast=contrast)

# Obtener los resultados del análisis diferencial
top_genes <- topTags(res, n=10)

# Imprimir los microARNs más diferencialmente expresados
top_genes$table
