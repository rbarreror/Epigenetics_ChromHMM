#!/bin/bash

# Autor: Rafael Barrero Rodríguez
# Fecha: 2020-03-12
# Descripción: Con este script haremos los distintos cálculos que se precisan para esta parte,
# y que consiste en analizar el solapamiento entre las regiones hiper e hipo metiladas y 
# y los segmentos asociados al estado 6.

SEGMENTS="segments_collapsed.bed"
HYPO="C001UYA3bs.hypo_meth.bs_call.20130415.bed"
HYPER="C001UYA3bs.hyper_meth.bs_call.20130415.bed"

LENGTH_SEG=$(awk -F"\t" 'BEGIN{suma=0};{suma+=$3-$2};END{print suma}' $SEGMENTS)
echo -e "\nNumero de bases en el estado 6: $LENGTH_SEG\n"

HYPO_BED=$(cut -f1-3 $HYPO | bedtools sort | bedtools merge)
HYPER_BED=$(cut -f1-3 $HYPER | bedtools sort | bedtools merge)


# Primero veremos qué porcentaje de nuestros segmentos podemos "capturar"
# si consideramos las regiones hypermetiladas

INTERSECT_HYPER=$(bedtools intersect -a $SEGMENTS -b <(echo -e "$HYPER_BED") | bedtools sort | bedtools merge)
PERC_HYPER=$(awk -v total=$LENGTH_SEG -F"\t" 'BEGIN{suma=0}; {suma+=$3-$2}; END{print 100*suma/total}' <(echo -e "$INTERSECT_HYPER"))

echo -e "Porcentaje de segmentos en regiones Hiper-metiladas: $PERC_HYPER"


# A continuación vemos qué porcentaje de nuestros segmentos podemos "capturar"
# si consideramos las regiones hypometiladas

INTERSECT_HYPO=$(bedtools intersect -a $SEGMENTS -b <(echo -e "$HYPO_BED") | bedtools sort | bedtools merge)
PERC_HYPO=$(awk -v total=$LENGTH_SEG -F"\t" 'BEGIN{suma=0}; {suma+=$3-$2}; END{print 100*suma/total}' <(echo -e "$INTERSECT_HYPO"))

echo -e "Porcentaje de segmentos en regiones Hipo-metiladas: $PERC_HYPO"

# Hacemos ahora el análisis inverso. Si consideramos nuestros segmentos, qué
# porcentaje de las regiones hiper e hipo metiladas capturamos

LENGTH_HYPER=$(awk -F"\t" 'BEGIN{suma=0};{suma+=$3-$2};END{print suma}' <(echo -e "$HYPER_BED"))
echo -e "\nNumero de bases en regiones Hypermetiladas: $LENGTH_HYPER"
PERC_HYPER_INV=$(awk -v total=$LENGTH_HYPER -F"\t" 'BEGIN{suma=0}; {suma+=$3-$2}; END{print 100*suma/total}' <(echo -e "$INTERSECT_HYPER"))
echo -e "Porcentaje de regiones Hiper-metiladas en nuestro estado: $PERC_HYPER_INV"



LENGTH_HYPO=$(awk -F"\t" 'BEGIN{suma=0};{suma+=$3-$2};END{print suma}' <(echo -e "$HYPO_BED"))
echo -e "\nNumero de bases en regiones Hypometiladas: $LENGTH_HYPO"
PERC_HYPO_INV=$(awk -v total=$LENGTH_HYPO -F"\t" 'BEGIN{suma=0}; {suma+=$3-$2}; END{print 100*suma/total}' <(echo -e "$INTERSECT_HYPO"))
echo -e "Porcentaje de regiones Hipo-metiladas en nuestro estado: $PERC_HYPO_INV"


# A continuación calculamos el % de segmentos en CpG que capturamos considerando
# regiones hipometiladas P(Hypo | E6, CpG)

SEGMENTS_CpG="intersect_cpg_ext.bed"

INTERSECT_HYPO_S_CpG=$(bedtools intersect -a <(bedtools sort -i $SEGMENTS_CpG | bedtools merge) -b <(echo -e "$HYPO_BED") | bedtools sort | bedtools merge)

LENGTH_S_CpG=$(awk -F"\t" 'BEGIN{suma=0};{suma+=$3-$2};END{print suma}' $SEGMENTS_CpG)

PERC_HYPO_S_CpG=$(awk -v total=$LENGTH_S_CpG -F"\t" 'BEGIN{suma=0}; {suma+=$3-$2}; END{print 100*suma/total}' <(echo -e "$INTERSECT_HYPO_S_CpG"))
echo -e "\nPorcentaje de segmentos en CpG que se encuentran en regiones Hypometiladas: $PERC_HYPO_S_CpG"

# Hacemos lo mismo pero con regiones hiper metiladas. Es decir, vemos qué porcentaje de los
# segmentos E6 que están en CpG se encuentran en regiones hiper-metiladas

INTERSECT_HYPER_S_CpG=$(bedtools intersect -a <(bedtools sort -i $SEGMENTS_CpG | bedtools merge) -b <(echo -e "$HYPER_BED") | bedtools sort | bedtools merge)

PERC_HYPER_S_CpG=$(awk -v total=$LENGTH_S_CpG -F"\t" 'BEGIN{suma=0}; {suma+=$3-$2}; END{print 100*suma/total}' <(echo -e "$INTERSECT_HYPER_S_CpG"))
echo -e "Porcentaje de segmentos en CpG que se encuentran en regiones Hypermetiladas: $PERC_HYPER_S_CpG"






















