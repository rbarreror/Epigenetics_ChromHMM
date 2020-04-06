#!/bin/bash

# Autor: Rafael Barrero Rodríguez
# Fecha: 2020-03-10
# Descripción: Cálculo a partir de las intersecciones del porcentaje de segmentos
# presente en la diferentes regiones funcionales del genoma.


FOLDER=$1
SEGMENTS=$2

BASES_STATE=$(awk -F "\t" 'BEGIN {suma = 0}; {suma += $3-$2}; END {print suma}' $SEGMENTS)

for FILE in $(ls $FOLDER/*.bed); do
    PERC=$(awk -v state=$BASES_STATE -F "\t" 'BEGIN {suma = 0}; {suma += $3-$2}; END {print(100*suma/state)}' $FILE)
    echo "$FILE(%): $PERC"
done
