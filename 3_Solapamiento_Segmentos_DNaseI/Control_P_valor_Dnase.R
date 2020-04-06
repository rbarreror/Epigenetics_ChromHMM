"""
Autor: Rafael Barrero Rodríguez
Fecha: 2020-03-11
Descripción: Análisis del p-valor en el fichero narrowPeaks. De este
modo podremos determinar si el fichero ha sido filtrado, procesado y normalizado.
"""


# Con este script comprobamos la calidad de los ficheros de DNasa

root_folder = paste0("/home/rafael/Master_UAM/Transcriptomica_RegulacionGenomica_Epigenomica/",
                     "3_Regulacion_Genomica_Epigenomica/Trabajo_Polycomb/")


setwd(paste0(root_folder, "DnaseI/blood_cells"))
ficheros_peaks <- list.files()
ficheros_peaks <- ficheros_peaks[grep(".narrowPeak", ficheros_peaks)]

narrow_peak_df <- read.table(ficheros_peaks[7], sep = "\t", header = FALSE)

# Representamos en un histograma la distribución del p-valor
p_values <- 10^(-narrow_peak_df[, 8])
hist(p_values, breaks = 20, col = "lightblue", xlab = "P-value",
     ylab = "Frequency", main = "P-value Distribution")


# Parece que ya han filtrado con un -log(pvalor) = 1.3, que equivale
# a un p-valor de 0.05 aproximadamente. 

# Comprobamos distribución de los scores, para ver si filtrar por aquí
# Representamos en un plot el p-valor frente al score
score <- narrow_peak_df[, 5]
hist(score, breaks = 20, col = "lightblue", xlab = "Score",
     ylab = "Frequency", main = "Score Distribution", xaxt = "n")
axis(side=1, at=seq(550, 1000, 50), labels = FALSE)
text(x=seq(600, 1000, 100),  -3000, 
     labels = seq(600, 1000, 100), pos = 1, xpd = TRUE)
plot(-log10(p_values), score, xlab = "-log(p-value)", ylab = "Score",
     main = "-log(p-value) vs Score", col = "blue", pch = 16, cex = 0.3)

# Comprobamos la distribución de la señal de los picos
enrichment <- narrow_peak_df[, 7]
hist(enrichment, breaks = 30, col = "lightblue", xlab = "Enrichment",
     ylab = "Frequency", main = "Enrichment Distribution")

plot(-log10(p_values), enrichment, xlab = "-log(p-value)", ylab = "Enrichment",
     main = "-log(p-value) vs Enrichment", col = "blue", pch = 16, cex = 0.3)

