#!/bin/bash

# Autor: Rafael Barrero Rodríguez
# Fecha: 2020-03-13
# Descripción: Con este script calculamos cómo varía el porcentaje de nucleótidos de nuestro estado
# conforme aumentamos la amplitud del intervalo en torno al TSS


# Con este comando generamos el BED con la posición del TSS
#cat TSS.mart | awk -F "\t" 'BEGIN {OFS="\t"}; NR!=1 {print("chr"$1, $2-1, $2)}' | bedtools sort | bedtools merge > TSS_mart.bed

# Fichero BED con los TSS
TSS_init=$1

#Fichero BED con los segmentos
SEGMENTS=$2

# ¿Analiza enhancers o CpG?
ENHANCER=1

if [ $ENHANCER -eq 1 ]; then
   TOTAL_B=7311600
else 
   TOTAL_B=$(cat $SEGMENTS | awk -F "\t" 'BEGIN {suma=0}; {suma += $3-$2}; END {print suma}')
fi
echo "Bases totales de nuestro estado: $TOTAL_B"

# Generamos fichero con porcentajes ampliando en una sola dirección
echo "" > percentage_unidir.tsv

# Generamos fichero con porcentajes ampliando en ambas direcciones
echo "" > percentage_bidir.tsv


for i in $(seq 0 20 2500); do
    echo "$i"
    # TSS.bed ampliando ANTES del gen
    TSS_MOD1=$(cat $TSS_init | awk -v interval=$i -F "\t" 'BEGIN {OFS="\t"}; {if(($2-interval)<0) print($1, 0, $3); else print($1, $2-interval, $3);}' | bedtools sort | bedtools merge)
    INTERSECT1=$(bedtools intersect -a $SEGMENTS -b <(echo -e "$TSS_MOD1") | bedtools sort | bedtools merge)
    PERC1=$(cat <(echo -e "$INTERSECT1") | awk -v total=$TOTAL_B -F "\t" 'BEGIN {suma=0}; {suma += $3-$2}; END {print(100*suma/total)}')
    echo -e "-$i\t$PERC1" >> percentage_unidir.tsv
    
    # TSS.bed ampliando DESPUES del gen
    TSS_MOD2=$(cat $TSS_init | awk -v interval=$i -F "\t" 'BEGIN {OFS="\t"}; {print($1, $2, $3+interval)}' | bedtools sort | bedtools merge)
    INTERSECT2=$(bedtools intersect -a $SEGMENTS -b <(echo -e "$TSS_MOD2") | bedtools sort | bedtools merge)
    PERC2=$(cat <(echo -e "$INTERSECT2") | awk -v total=$TOTAL_B -F "\t" 'BEGIN {suma=0}; {suma += $3-$2}; END {print(100*suma/total)}')
    echo -e "$i\t$PERC2" >> percentage_unidir.tsv

    # TSS.bed ampliando en ambas direcciones
    TSS_MOD3=$(cat $TSS_init | awk -v interval=$i -F "\t" 'BEGIN {OFS="\t"}; {if(($2-interval)<0) print($1, 0, $3+interval); else print($1, $2-interval, $3+interval);}' | bedtools sort | bedtools merge)
    INTERSECT3=$(bedtools intersect -a $SEGMENTS -b <(echo -e "$TSS_MOD3") | bedtools sort | bedtools merge)
    PERC3=$(cat <(echo -e "$INTERSECT3") | awk -v total=$TOTAL_B -F "\t" 'BEGIN {suma=0}; {suma += $3-$2}; END {print(100*suma/total)}')
    echo -e "$i\t$PERC3" >> percentage_bidir.tsv
done
