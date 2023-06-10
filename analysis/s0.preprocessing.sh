#!/bin/bash

#PBS -J 1-22
#PBS -k oed
#PBS -j oe
#PBS -l nodes=1:ppn=1 -l mem=1G
#PBS -V
#PBS -S /bin/bash
#PBS -q workq
#PBS -o output/anno_vcfs/logs/parse_mbiobank_log.jobarray.o 
#PBS -e output/anno_vcfs/logs/parse_mbiobank_log.jobarray.e

cd $PBS_O_WORKDIR
source ~/miniconda3/bin/activate
conda activate vep
mkdir -p output/anno_vcfs

n_count=22
if [[ $(($PBS_ARRAY_INDEX % 22)) -ne 0 ]];then
    n_count=$(($PBS_ARRAY_INDEX % 22))
fi

invcf=$(grep mbiobank all_vcfs.list |head -n $n_count | tail -n1)
vcffn=$(basename $invcf)
annovcf=output/anno_vcfs/raw_anno/${vcffn%.vcf.gz}.anno.vcf.gz
tmpvcf=output/anno_vcfs/raw_anno/${vcffn%.vcf.gz}.anno.tmp.vcf.gz

echo ">PBS_ARRAY_INDEX:$PBS_ARRAY_INDEX"
echo ">invcf:$invcf"
echo ">vcffn:$vcffn"
echo ">annovcf:$annovcf"
echo ">tmpvcf:$tmpvcf"

[[ -f $tmpvcf ]] && echo "$annovcf not finished" && exit 

POP=CHN

mkdir -p output/anno_vcfs/mbiobank/$POP
outvcf=output/anno_vcfs/mbiobank/$POP/${vcffn%.vcf.gz}.anno.vcf.gz
passvcf=output/pass_vcfs/mbiobank/$POP/${vcffn%.vcf.gz}.anno.$POP.vcf.gz
echo "Generating $outvcf..."

set -x
/ehpcdata/analysis/linzhe/lowDepth/resourse/software/bcftools-1.16/bcftools view -e "AC == 0 | AF > 0.05" -Oz -o $passvcf $outvcf
/ehpcdata/analysis/linzhe/lowDepth/resourse/software/bcftools-1.16/bcftools index $passvcf
/ehpcdata/analysis/linzhe/lowDepth/resourse/software/bcftools-1.16/bcftools stats $passvcf > $passvcf.stats
set +x
