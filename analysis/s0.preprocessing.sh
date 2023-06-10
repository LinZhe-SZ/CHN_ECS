#!/bin/bash
set -eof pipefail

# ChinaMap
vcf=ChinaMAP.phase1.vcf.gz
outdir=output/pass_vcfs/ChinaMap/
mkdir -p outdir

for i in {1..22};do
    passvcf=$outdir/ChinaMap.chr$i.pass.vcf.gz
    bcftools index -f $vcf
    bcftools norm -m - $vcf chr$i | bcftools view -e "AC==0 | AF > 0.05 | QUAL <= 100 | AN < 10229" -Oz -o $passvcf 
    bcftools index $passvcf
    bcftools stats $passvcf > $passvcf.stats
done

# WBBC
vcflist=wbbc_vcfs.list # a list contains file path of WBBC VCFs, one file per line
outdir=output/pass_vcfs/WBBC/
mkdir -p outdir
cat $vcflist | while read -r vcf;do
    passvcf=$outdir/WBBC.chr$i.pass.vcf.gz
    ln $vcf $passvcf # all variants were considered high quality
    bcftools index $passvcf
    bcftools stats $passvcf > $passvcf.stats
done

# gnomAD v3
vcflist=gnomAD_vcfs.list # a list contains file path of gnomAD VCFs, one file per line

POPS=("EAS" "EUR" "AFR" "AMR" "SAS")
for pop in "${POPS[@]}";do # extract variants of each sub-population respectively
    outdir=output/pass_vcfs/gnomAD_${pop}/
    mkdir -p outdir
    cat $vcflist | while read -r vcf;do
        pn=$(echo $POP | tr '[:upper:]' '[:lower:]')
        if [[ $pn == "eur" ]];then
            pn='nfe'
        fi
        passvcf=$outdir/gnomAD_${pop}.chr$i.pass.vcf.gz
        bcftools norm -m - $vcf | bcftools view -e "AF_non_cancer_${pn} == 0 | AF_non_cancer_${pn} > 0.05" -i "FILTER=="PASS"' -Oz -o $passvcf
        bcftools index $passvcf
        bcftools stats $passvcf > $passvcf.stats
    done
done
