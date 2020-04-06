"""
Autor: Rafael Barrero Rodríguez
Fecha: 2020-03-15
Descripción: Con este script realizamos un doble filtrado de los segmentos, por intersección
entre ambas réplicas, y por probabilidad posterior. Basándonos en análsis del script
defining_threshold.R, decidimos filtrar usando 0.75 como umbral.
"""
root_folder = paste0("/home/rafael/Master_UAM/Transcriptomica_RegulacionGenomica_Epigenomica/",
                     "3_Regulacion_Genomica_Epigenomica/Trabajo_Polycomb/")

estado <- 6

# Ejercicio 1 - Encontrar segmentos con el Estado 6 en ambos replicados

# Cargar ficheros con los estados

setwd(paste0(root_folder, "RESULTS/Modelo_11_estados/STATEBYLINE"))
ficheros_state_by_line <- list.files()
ficheros_state_by_line <- ficheros_state_by_line[order(ficheros_state_by_line)]
chrs <- paste0(rep("_chr", 22), as.character(seq(1, 22)), "_")
chrs <- c(chrs, c("_chrM_", "_chrX_"))

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


# Filtrar por probabilidad posterior (si se quiere)
setwd(paste0(root_folder, "RESULTS/Modelo_11_estados/POSTERIOR"))
ficheros_posterior_prob <- list.files()
ficheros_posterior_prob <- ficheros_posterior_prob[order(ficheros_posterior_prob)]

threshold <- 0.75

# Esta función recibe el número de un cromosoma ("_chrN") y lee para cada monocito
# el fichero con la probabilidad posterior, quedándose con la columna del estado
# que nos interesa (el 6). Devuelve un booleano con TRUE en los segmentos que tienen
# una probabilidad superior al threshold en ambos monocitos.
get_segments_above_threshold <- function(chr){
  ficheros_chr <- ficheros_posterior_prob[grep(chr, ficheros_posterior_prob)]
  fichero_mon1 <- ficheros_chr[grep("Monocyte1", ficheros_chr)]
  fichero_mon2 <- ficheros_chr[grep("Monocyte2", ficheros_chr)]
  posterior_prob_mon1 <- read.table(fichero_mon1, sep = "\t", header = TRUE, skip = 1)[, estado]
  posterior_prob_mon2 <- read.table(fichero_mon2, sep = "\t", header = TRUE, skip = 1)[, estado]
  posterior_prob_mon1_above <- posterior_prob_mon1 > threshold
  posterior_prob_mon2_above <- posterior_prob_mon2 > threshold
  return(posterior_prob_mon1_above & posterior_prob_mon2_above)
}

bool_common_post_prob <- sapply(chrs, FUN = get_segments_above_threshold)


# A continuación nos quedamos con los segmentos que tienen el estado en ambos monocitos
# y que superan el threshold. Esta operación es redundante, pues si el threshold es de 0.8
# el resultado será igual al de bool_common_post_prob.

get_segments_with_common_state_and_post_prob <- function(chr){
  return(bool_common_state[[chr]] & bool_common_post_prob[[chr]])
}

bool_common_state_and_post_prob <- sapply(chrs, FUN = get_segments_with_common_state_and_post_prob)

# Para hacer el fichero bed necesitamos las posiciones de inicio y fin de cada segmento.
# A continuación obtenemos vectores con el inicio de cada segmento por cada cromosoma

# Avoid exponential notation
options(scipen=999)

get_segment_start <- function(chr){
  n_segments <- length(bool_common_state[[chr]])
  vec <- seq(0, 200*(n_segments - 1), 200)
  return(vec)
}

segment_start <- sapply(chrs, FUN = get_segment_start)

# A continuación obtenemos vectores con el final de cada segmento por cada cromosoma

get_segment_end <- function(chr){
  n_segments <- length(bool_common_state[[chr]])
  vec <- seq(200, 200*(n_segments), 200)
  return(vec)
}

segment_end <- sapply(chrs, FUN = get_segment_end)

# Construimos data frame a partir del cual generamos el bed

pre_bed <- data.frame(chr = c(), start = c(), end = c())

for(chr in chrs){
  state_of_interest <- bool_common_state_and_post_prob[[chr]]
  state_of_interest_index <- which(state_of_interest)
  
  start_col <- segment_start[[chr]][state_of_interest_index]
  end_col <- segment_end[[chr]][state_of_interest_index]
  
  chr_name <- stringr::str_remove_all(chr, "_")
  chr_col <- rep(chr_name, length(start_col))
  
  pre_bed <- rbind.data.frame(pre_bed, data.frame(chr = chr_col, start = start_col,
                                                  end = end_col))
}

# El pre_bed lo escribimos en un fichero
setwd(root_folder)
write.table(pre_bed, file = paste0(root_folder, "output_R/segments2.bed"), quote = FALSE,
            col.names = FALSE, row.names = FALSE, sep = "\t")
