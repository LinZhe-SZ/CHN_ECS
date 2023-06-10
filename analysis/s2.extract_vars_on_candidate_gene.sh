#!/bin/bash
set -eof pipefail

for vcf in output/anno_vcfs/*/*vcf.gz;do
    outvcf=$(sed 's/anno_vcfs/candidate_vars/' <<< $vcf)
    mkdir -p $(dirname $outvcf)

    python scripts/parse_vcf_by_gene.py -i $vcf -g data/candidate_gene_excludeXY.list -o $outvcf > $outvcf.parse.log
    bcftools index -f $outvcf
    bcftools stats $outvcf > $outvcf.stats

done
