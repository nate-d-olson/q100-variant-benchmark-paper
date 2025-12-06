#!/usr/bin/bash

VCF=$1
BED=$2

 bcftools query -HH \
            -f '%CHROM\t%POS0\t%END\t%SVTYPE\t%VKX\t%SVLEN\t%TRF\t%TRFstart\t%TRFend\t%RM_clsfam\n' \
            -R ${BED} ${VCF}
