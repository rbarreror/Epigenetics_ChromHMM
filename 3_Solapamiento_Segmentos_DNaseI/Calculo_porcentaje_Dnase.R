"""
Autor: Rafael Barrero Rodríguez
Fecha: 2020-03-11
Descripción: Calculamos los porcentajes de solapamientos entre
las regiones desplegadas según el fichero narrowPeaks y los
nucleótidos presentes en estado bivalente
"""

# Cálculo del porcentaje de coincidencia con DNaseI
# Lo calculamos como Bases de Estado en DNaseI / Bases de Estado Total

root_folder = paste0("/home/rafael/Master_UAM/Transcriptomica_RegulacionGenomica_Epigenomica/",
                     "3_Regulacion_Genomica_Epigenomica/Trabajo_Polycomb/")

col_names <- c("chr", "start", "end")


# Abrir fichero con los fragmentos de DNaseI
dnaseI_interval <- read.table(paste0(root_folder, 
                                     "DnaseI/wgEncodeOpenChromDnaseMonocd14Pk.narrowPeak"))[, 1:3]
colnames(dnaseI_interval) <- col_names

# Abrir fichero con los segmentos colapsados de nuestro estado
state_bed_interval <- read.table(paste0(root_folder, "DnaseI/segments_collapsed.bed"))
colnames(state_bed_interval) <- col_names

# Abrir fichero con las intersecciones obtenidas con bedtools
intersect_bed_interval <- read.table(paste0(root_folder, 
                                            "DnaseI/intersect/intersect_segments_14Pk.bed"))
colnames(intersect_bed_interval) <- col_names


# Calculamos número de bases donde hay intersección

get_interval_dif <- function(intersect_interval){
  start_interval <- as.integer(intersect_interval[2])
  end_interval <- as.integer(intersect_interval[3])
  return(end_interval - start_interval)
}

base_pair_intersect <- sum(apply(intersect_bed_interval, MARGIN = 1, get_interval_dif))
base_pair_state <- sum(apply(state_bed_interval, MARGIN = 1, get_interval_dif))
base_pair_dnaseI <- sum(apply(dnaseI_interval, MARGIN = 1, get_interval_dif))

# Porcentaje de bases de nuestro estado que aparecen en fragmentos de DNaseI
100 * base_pair_intersect/base_pair_state


# Porcentaje de bases de fragmento que aparecen en nuestro estado
100 * base_pair_intersect/base_pair_dnaseI
