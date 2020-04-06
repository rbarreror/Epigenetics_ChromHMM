"""
Autor: Rafael Barrero Rodríguez
Fecha: 2020-03-05
Descripción: Analizamos cómo se distribuye la probabilidad posterior de pertenencia
al estado 6 en cada uno de los segmentos P(E6|Segmento). Asimismo, vemos el
porcentaje de segmentos eliminados en función del umbral establecido
"""


# Script para obtener gráfica con probabilidades posteriores y porcentaje de excluidos

root_folder = paste0("/home/rafael/Master_UAM/Transcriptomica_RegulacionGenomica_Epigenomica/",
                     "3_Regulacion_Genomica_Epigenomica/Trabajo_Polycomb/")

estado <- 6

setwd(paste0(root_folder, "RESULTS/Modelo_11_estados/POSTERIOR"))
ficheros_posterior_prob <- list.files()
ficheros_posterior_prob <- ficheros_posterior_prob[order(ficheros_posterior_prob)]

threshold <- 0.7

chrs <- paste0(rep("_chr", 22), as.character(seq(1, 22)), "_")
chrs <- c(chrs, c("_chrM_", "_chrX_"))

# Cargar ficheros con los estados

setwd(paste0(root_folder, "RESULTS/Modelo_11_estados/STATEBYLINE"))
ficheros_state_by_line <- list.files()
ficheros_state_by_line <- ficheros_state_by_line[order(ficheros_state_by_line)]

# La función recibe el número de cromosoma ("_chrN_"), y lee los estados
# para cada monocito. Devuelve un booleano donde coincidan que los estados
# son igual a 6 en ambos monocitos
get_common_state_line <- function(chr){
  ficheros_chr <- ficheros_state_by_line[grep(chr, ficheros_state_by_line)]
  fichero_mon1 <- ficheros_chr[grep("Monocyte1", ficheros_chr)]
  fichero_mon2 <- ficheros_chr[grep("Monocyte2", ficheros_chr)]
  states_mon1 <- as.integer(readLines(fichero_mon1)[-c(1,2)])
  states_mon1_bool <- states_mon1 == estado
  states_mon2 <- as.integer(readLines(fichero_mon2)[-c(1,2)])
  states_mon2_bool <- states_mon2 == estado
  return(states_mon1_bool & states_mon2_bool)
}

# Contiene una lista de vectores booleanos. Cada vector se corresponde a los segmentos
# de un cromosoma y tiene un TRUE en el segmento donde ambos monocitos tenían estado 6
bool_common_state <- sapply(chrs, FUN = get_common_state_line)



# La función toma las probabilidades posteriores de los segmentos donde aparece nuestro estado
# en ambas réplicas. Para cada segmento tiene dos probabilidades posteriores, de modo que se
# queda con la probabilidad mínima. Devuelve un vector donde para cada segmento de nuestro estado
# contiene la probabilidad posterior mínima de entre ambas réplicas
setwd(paste0(root_folder, "RESULTS/Modelo_11_estados/POSTERIOR"))
ficheros_posterior_prob <- list.files()
ficheros_posterior_prob <- ficheros_posterior_prob[order(ficheros_posterior_prob)]

get_segments_post_prob <- function(chr){
  ficheros_chr <- ficheros_posterior_prob[grep(chr, ficheros_posterior_prob)]
  fichero_mon1 <- ficheros_chr[grep("Monocyte1", ficheros_chr)]
  fichero_mon2 <- ficheros_chr[grep("Monocyte2", ficheros_chr)]
  posterior_prob_mon1 <- read.table(fichero_mon1, sep = "\t", header = TRUE, skip = 1)[, estado]
  posterior_prob_mon1 <- posterior_prob_mon1[bool_common_state[[chr]]]
  posterior_prob_mon2 <- read.table(fichero_mon2, sep = "\t", header = TRUE, skip = 1)[, estado]
  posterior_prob_mon2 <- posterior_prob_mon2[bool_common_state[[chr]]]
  posterior_prob_joint <- apply(rbind(posterior_prob_mon1, posterior_prob_mon2), 2, min)
  return(posterior_prob_joint)
}

post_prob_chr <- sapply(chrs, FUN = get_segments_post_prob)

# Obtenemos un vector con las probabilidades posteriores de todos los cromosomas

all_posterior <- unlist(post_prob_chr)
hist(all_posterior, breaks = 25, col = "lightblue", xlab = "Posterior Probability",
     ylab = "Frequency", main = "Posterior Probability Distribution")

# Porcentaje con Posterior mayor que:

filter_ev <- sapply(seq(0, 1, 0.01), FUN = function(x) sum(all_posterior > x)/length(all_posterior))
plot(seq(0, 1, 0.01), 100-100*filter_ev, ylab = "Percentage of segments filtered", xlab = "Threshold",
     main = "% Segments filtered vs Threshold", type = "l", col = "blue", xaxt = "n", yaxt = "n", lwd = 1.5)
axis(side=1, at=seq(0, 1, 0.05), labels = FALSE)
text(x=seq(0, 1, 0.1),  -5, 
     labels = seq(0, 1, 0.1), pos = 1, xpd = TRUE)
axis(side=2, at=seq(0, 100, 10), labels = FALSE)
text(y=seq(0, 100, 10),  -0.05, 
     labels = seq(0, 100, 10), pos = 2, xpd = TRUE)

# Con un threshold del 0.75 excluimos un 10% de los segmentos 
sum(all_posterior > 0.75)/length(all_posterior)
