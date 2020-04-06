"""
Name: Rafael Barrero Rodríguez
Date: 2020-03-12
Description: En este script vamos a ver cómo se distribuyen los porcentajes
de metilación en los ficheros BED. El fichero indica para cada intervalo, el
porcentaje de veces que aparece metilado. Vamos a ver cuál esla distribución
de ese porcentaje en los distintos intervalos para asegurarnos de que no hay
que filtrar
"""

root_folder = paste0("/home/rafael/Master_UAM/Transcriptomica_RegulacionGenomica_Epigenomica/",
                     "3_Regulacion_Genomica_Epigenomica/Trabajo_Polycomb/metilacion/")

setwd(root_folder)
hyper_met_table <- read.table("C001UYA3bs.hyper_meth.bs_call.20130415.bed", 
                              header = FALSE, sep = "\t")
hypo_met_table <- read.table("C001UYA3bs.hypo_meth.bs_call.20130415.bed", 
                              header = FALSE, sep = "\t")


# Vemos distribucion de % en fichero con BEDs hipermetilados

perc_hyper <- hyper_met_table$V5
hist(perc_hyper, breaks = 20, col = "lightblue", xlab = "Methylation ratio",
     ylab = "Frequency", main = "Methylation ratio Distribution (Hyper)")

"""
Como vemos, todos los intervalos contenidos en el hyper BED presentan un ratio 
de metilación superior al 0.75. Por lo tanto, podemos considerar todos estos
intervalos como regiones hipermetiladas
"""


# Vemos distribucion de % en fichero con BEDs hipometilados

perc_hypo <- hypo_met_table$V5
hist(perc_hypo, breaks = 20, col = "lightblue", xlab = "Methylation ratio",
     ylab = "Frequency", main = "Methylation ratio Distribution (Hypo)")


"""
Como vemos, todos los intervalos contenidos en el hypo BED presentan un ratio 
de metilación inferior a 0.25. Por lo tanto, podemos considerar todos estos
intervalos como regiones hipometiladas
"""