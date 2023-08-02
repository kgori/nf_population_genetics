#!/bin/bash -ue
bcftools view --threads $3 -m 2 -M 2 $1 | \
  bcftools filter --threads $3 \
    -i '(((REF="A"|REF="G") & (ALT="C"|ALT="T")) | ((REF="C"|REF="T") & (ALT="A"|ALT="G")))' | \
    bcftools filter -e '(REF="A" & ALT="G") | (REF="G" & ALT="A") | (REF="C" & ALT="T") | (REF="T" & ALT="C")' \
      -Oz -o $2
