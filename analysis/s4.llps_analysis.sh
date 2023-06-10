#!/binb/bash
set -eof pipefail

outdir=ChinaMap_llps
mkdir -p $outdir

cat output/GCR/ChinaMap/*.pathogenic.tsv | awk 'BEGIN{OFS="\t"}NR>1{print $1,$2-1,$2}' > $outdir/pathogenic.var.bed

for i in {1..22};do
    bcftools view -T $outdir/pathogenic.var.bed -Oz -o $outdir/ChinaMap.chr$i.llps-pathogenic.vcf.gz
    bcftools +split-vep -c Uniprot_acc,Consequence,HGVSp_VEP,VEP_canonical,SYMBOL,Ensembl_proteinid,Amino_acids,Codons, -f "%CHROM\t%POS\t%REF\t%ALT\t%SYMBOL\t%Consequence\t%AF\t%Uniprot_acc\t%HGVSp_VEP\t%VEP_canonical\t%CLNSIG\t%CLNSIGCONF\n" $outdir/ChinaMap.chr$i.llps-pathogenic.vcf.gz >> $outdir/pathogenic.pro.list
done

python scripts/llps_annotation.py -i $outdir/pathogenic.pro.list -l data/1-s2.0-S1534580722004506-mmc2.xls -o $outdir/pathogenic.pro.anno.list

