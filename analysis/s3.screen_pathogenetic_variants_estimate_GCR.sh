#!/bin/bash

for vcf in output/candidate_vars/*/*vcf.gz;do
    vcffn=$(basename $vcf)
    db=$(basename $(dirname $vcf))
    mkdir -p output/$db

    outtable=output/$db/${vcffn%.vcf.gz}.table
    outgcr=output/$db/${vcffn%.vcf.gz}.gcr

    if [[ $db == "ChinaMap" || $db == "WBBC" ]];then
        addparams="--an-name AN --ac-name AC --af-name AF"
    else
        pop=$(cut -f2 -d'_' <<< $db | tr '[:upper:]' '[:lower:]')
        [[ $pop == "eur" ]] && pop="nfe"
        addparams="--an-name AN_non_cancer_$pop --ac-name AC_non_cancer_$pop --af-name AF_non_cancer_$pop"
    fi

    python scripts/screen_pathogenetic_variant.py $addparams -f -v -i $vcf -p $outtable -g $outgcr
done