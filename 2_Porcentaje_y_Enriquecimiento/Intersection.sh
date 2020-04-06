#!/bin/bash

# Autor: Rafael Barrero Rodríguez
# Fecha: 2020-03-10
# Descripción: En este script hacemos las intersecciones entre nuestros segmentos y
# las distintas regiones. Lo guardaremos en la carpeta intersect

SEGMENTS=$1
SEGMENTS_NAME=$(basename -s .bed $SEGMENTS)

bedtools intersect -a $SEGMENTS -b <(bedtools sort -i  genes.bed) | bedtools sort | bedtools merge > intersect/intersect_${SEGMENTS_NAME}_genes.bed

bedtools intersect -a $SEGMENTS -b <(bedtools sort -i TSS.bed) | bedtools sort | bedtools merge > intersect/intersect_${SEGMENTS_NAME}_TSS.bed

bedtools intersect -a $SEGMENTS -b <(bedtools sort -i coding_region.bed) | bedtools sort | bedtools merge > intersect/intersect_${SEGMENTS_NAME}_coding.bed

bedtools intersect -a $SEGMENTS -b <(bedtools sort -i enhancers_HACER_VISTA.bed) | bedtools sort | bedtools merge > intersect/intersect_${SEGMENTS_NAME}_enhancers.bed

bedtools intersect -a $SEGMENTS -b <(bedtools sort -i cpgIslandExt.hg19.bed) | bedtools sort | bedtools merge > intersect/intersect_${SEGMENTS_NAME}_CpG.bed

bedtools intersect -a $SEGMENTS -b <(bedtools sort -i exon.bed) | bedtools sort | bedtools merge > intersect/intersect_${SEGMENTS_NAME}_exon.bed

bedtools intersect -a $SEGMENTS -b <(bedtools sort -i TSS_200bp.bed) | bedtools sort | bedtools merge > intersect/intersect_${SEGMENTS_NAME}_TSS_200.bed

bedtools intersect -a $SEGMENTS -b <(bedtools sort -i TSS_2000bp.bed) | bedtools sort | bedtools merge > intersect/intersect_${SEGMENTS_NAME}_TSS_2000.bed

