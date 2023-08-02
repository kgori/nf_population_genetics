#!/bin/bash -ue
bcftools query --format '%CHROM\t%POS\t%REF\t%ALT\n' $1 | awk 'BEGIN {OFS="\t"} {print $1,$2-1,$2,NR"-"$3"_"$4}' | bgzip > id_annotation.bed.gz
