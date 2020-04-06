"""
Autor: Rafael Barrero Rodríguez
Fecha: 2020-03-13
Descripción: Con este script representamos cómo varía el porcentaje de nuestro estado a medida
que ampliamos el intervalo de TSS. Es una representación de los datos obtenidos con el script
Distribucion_TSS.sh
"""


root_folder = paste0("/home/rafael/Master_UAM/Transcriptomica_RegulacionGenomica_Epigenomica/",
                     "3_Regulacion_Genomica_Epigenomica/Trabajo_Polycomb/anotation/",
                     "distribucionBasica/")
setwd(root_folder)


tabla_unidir <- read.table("TSS_ensembl/resultado_0_2500_10/percentage_unidir.tsv", 
                           header = FALSE, sep = "\t")
tabla_bidir <- read.table("TSS_ensembl/resultado_0_2500_10/percentage_bidir.tsv", 
                          header = FALSE, sep = "\t")

tss_variation <- tabla_bidir[, 1]
perc_pre <- tabla_unidir[tabla_unidir[, 1] <= 0, 2][-1]
perc_post <- tabla_unidir[tabla_unidir[, 1] >= 0, 2][-1]
perc_bi <- tabla_bidir[, 2]


# Representamos cómo varía el porcentaje de nucleótidos de nuestro estado
# a medida que ampliamos la ventana en torno al TSS. En rojo ampliamos la
# ventana hacia la zona reguladora (antes de la transcripción). En azul la
# ampliamos hacia la zona codificante. En verde la ampliamos en ambas direcciones
# El comportamiento es similar al ampliar el intervalo a una y otra dirección
# del TSS. La suma de los porcentajes de estos dos intervalos es mayor que 
# el porcentaje obtenido en el intervalo bidireccional debido al solapamiento
# entre regiones flanqueantes de distintos TSS.

# El porcentaje de saturación en los casos unidireccionales es del 50%, mientras
# que en el bidireccional es del 76%. Las tres curvas comienzan su saturación a 
# partir de los 1000 nucleótidos (de distancia del TSS) aproximadamente.

# En cualquier caso, estas gráficas revelan que nuestro estado se encuentra
# en un gran porcentaje en las zonas flanqueantes al inicio de la transcripción (+-1000 pb),
# próximas a las zonas reguladoras y promotoras.

plot(tss_variation, perc_pre, type="l", col="red", xlab = "Pair Bases from TSS",
     ylab = "Percentage of State 6", main = "Percentage of State 6 vs TSS interval", lwd = 1.5,
     ylim = c(0, 100), xaxt = "n")
lines(tss_variation, perc_post, col="blue")
lines(tss_variation, perc_bi, col="darkgreen")
axis(side=1, at=seq(0, 2500, 500), labels = TRUE)
legend(-80, 100, legend=c("Before TSS", "After TSS", "Both sides"),
       col=c("red", "blue", "darkgreen"), lty=1, cex=0.8, box.lty = 0)


#################################################
# Comparamos con TES
#################################################
tabla_unidir_tes <- read.table("TES_ensembl/resultado_0_2500_30/percentage_unidir.tsv", 
                           header = FALSE, sep = "\t")
tabla_bidir_tes <- read.table("TES_ensembl/resultado_0_2500_30/percentage_bidir.tsv", 
                          header = FALSE, sep = "\t")

tss_variation_tes <- tabla_bidir_tes[, 1]
perc_pre_tes <- tabla_unidir_tes[tabla_unidir_tes[, 1] <= 0, 2][-1]
perc_post_tes <- tabla_unidir_tes[tabla_unidir_tes[, 1] >= 0, 2][-1]
perc_bi_tes <- tabla_bidir_tes[, 2]


# Con este gráfico comparamos cómo varía la presencia de nuestro estado al
# abrir ventanas en TSS y en TES. Vemos que las ventanas abiertas en torno
# a TSS capturan un porcentaje mayor de nuestro estado que las abiertas en
# torno a TES.
plot(tss_variation, perc_bi, type="l", col="darkgreen", xlab = "Pair Bases from TSS",
     ylab = "Percentage of State 6", main = "Percentage of State 6 vs TSS/TES interval", lwd = 1.5,
     ylim = c(0, 100), xaxt = "n")
lines(tss_variation_tes, perc_bi_tes, col="darkblue", lwd = 1.5)
axis(side=1, at=seq(0, 2500, 500), labels = TRUE)
legend(-80, 100, legend=c("TSS", "TES"),
       col=c("darkgreen", "darkblue"), lty=1, cex=0.8, box.lty = 0)

# Podemos comparar Predecesor de TSS con Posterior de TES
plot(tss_variation, perc_pre, type="l", col="darkgreen", xlab = "Pair Bases from TSS",
     ylab = "% of State 6", main = "Percentage of State 6 vs TSS/TES interval", lwd = 1.5,
     ylim = c(0, 100), xaxt = "n")
lines(tss_variation_tes, perc_post_tes, col="darkblue", lwd = 1.5)
axis(side=1, at=seq(0, 20000, 2000), labels = TRUE)
legend(-80, 100, legend=c("Before TSS", "After TES"),
       col=c("darkgreen", "darkblue"), lty=1, cex=0.8, box.lty = 0)


#########################################################
# Repetimos los anterior pero con ventanas de hasta 20000
#########################################################

tabla_unidir <- read.table("TSS_ensembl/resultado_0_20000_400/percentage_unidir.tsv", 
                           header = FALSE, sep = "\t")
tabla_bidir <- read.table("TSS_ensembl/resultado_0_20000_400/percentage_bidir.tsv", 
                          header = FALSE, sep = "\t")

tss_variation <- tabla_bidir[, 1]
perc_pre <- tabla_unidir[tabla_unidir[, 1] <= 0, 2][-1]
perc_post <- tabla_unidir[tabla_unidir[, 1] >= 0, 2][-1]
perc_bi <- tabla_bidir[, 2]

plot(tss_variation, perc_pre, type="l", col="red", xlab = "Pair Bases from TSS",
     ylab = "% of State 6", main = "Percentage of State 6 vs TSS interval", lwd = 1.5,
     ylim = c(0, 100), xaxt = "n")
lines(tss_variation, perc_post, col="blue")
lines(tss_variation, perc_bi, col="darkgreen")
axis(side=1, at=seq(0, 20000, 2000), labels = TRUE)
legend(-600, 102, legend=c("Before TSS", "After TSS", "Both sides"),
       col=c("red", "blue", "darkgreen"), lty=1, cex=0.8, box.lty = 0)

# Comparamos con TES

tabla_unidir_tes <- read.table("TES_ensembl/resultado_0_20000_400/percentage_unidir.tsv", 
                               header = FALSE, sep = "\t")
tabla_bidir_tes <- read.table("TES_ensembl/resultado_0_20000_400/percentage_bidir.tsv", 
                              header = FALSE, sep = "\t")

tss_variation_tes <- tabla_bidir_tes[, 1]
perc_pre_tes <- tabla_unidir_tes[tabla_unidir_tes[, 1] <= 0, 2][-1]
perc_post_tes <- tabla_unidir_tes[tabla_unidir_tes[, 1] >= 0, 2][-1]
perc_bi_tes <- tabla_bidir_tes[, 2]


# Con este gráfico comparamos cómo varía la presencia de nuestro estado al
# abrir ventanas en TSS y en TES. Vemos que las ventanas abiertas en torno
# a TSS capturan un porcentaje mayor de nuestro estado que las abiertas en
# torno a TES.
plot(tss_variation, perc_bi, type="l", col="darkgreen", xlab = "Pair Bases from TSS",
     ylab = "% of State 6", main = "Percentage of State 6 vs TSS/TES interval", lwd = 1.5,
     ylim = c(0, 100), xaxt = "n")
lines(tss_variation_tes, perc_bi_tes, col="darkblue", lwd = 1.5)
axis(side=1, at=seq(0, 20000, 2000), labels = TRUE)
legend(-600, 100, legend=c("TSS", "TES"),
       col=c("darkgreen", "darkblue"), lty=1, cex=0.8, box.lty = 0)



# Podemos comparar Predecesor de TSS con Posterior de TES
plot(tss_variation, perc_pre, type="l", col="darkgreen", xlab = "Pair Bases from TSS",
     ylab = "% of State 6", main = "Percentage of State 6 vs TSS/TES interval", lwd = 1.5,
     ylim = c(0, 100), xaxt = "n")
lines(tss_variation_tes, perc_post_tes, col="darkblue", lwd = 1.5)
axis(side=1, at=seq(0, 20000, 2000), labels = TRUE)
legend(-80, 100, legend=c("Before TSS", "After TES"),
       col=c("darkgreen", "darkblue"), lty=1, cex=0.8, box.lty = 0)


###########################
# Análisis CpG
###########################

# Lo que hacemos a continuación es analizar los CpG que se encuentran
# en los alrededores del TSS. Tomamos en consideración los CpG que están
# en un intervalo alrededor del TSS, y qué porcentaje de segmentos
# logramos capturar. 

tabla_unidir <- read.table("TSS_ensembl/resultado_CpG_0_2500_20/percentage_unidir.tsv", 
                           header = FALSE, sep = "\t")
tabla_bidir <- read.table("TSS_ensembl/resultado_CpG_0_2500_20/percentage_bidir.tsv", 
                          header = FALSE, sep = "\t")

tss_variation <- tabla_bidir[, 1]
perc_pre <- tabla_unidir[tabla_unidir[, 1] <= 0, 2][-1]
perc_post <- tabla_unidir[tabla_unidir[, 1] >= 0, 2][-1]
perc_bi <- tabla_bidir[, 2]

plot(tss_variation, perc_pre, type="l", col="red", xlab = "Pair Bases from TSS",
     ylab = "% of State 6", main = "Percentage of State 6 in CpG", lwd = 1.5,
     ylim = c(0, 50), xaxt = "n")
lines(tss_variation, perc_post, col="blue")
lines(tss_variation, perc_bi, col="darkgreen")
axis(side=1, at=seq(0, 2500, 500), labels = TRUE)
legend(-80, 50, legend=c("Before TSS", "After TSS", "Both sides"),
       col=c("red", "blue", "darkgreen"), lty=1, cex=0.8, box.lty = 0)

# En este caso, vemos qué porcentaje de segmentos se encuentran en las CpG
# localizadas en un intervalo en torno al TSS. Conforme aumenta ese intervalo,
# vemos que un mayor % de segmentos se encuentra en las CpG dentro de ese
# intervalo. Tomando en consideración los CpG situados a una ventana de 2500
# bp en torno a los TSS, capturamos un 40% de los segmentos. Esto nos da una
# idea del papel regulador de nuestro estado.

