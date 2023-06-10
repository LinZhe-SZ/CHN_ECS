#!/bin/bash
set -eof pipefail

function anno () {
    VEPCACHEDIR=/path/to/VEP108_GRCh38_Cache/
    VEPPLUGINS=/path/to/ensembl-vep-108.1-0/
    DBNSFP=/path/to/dbNSFP/dbNSFP4.3a/dbNSFP4.3a.gz
    CLINVAR=/path/to/Clinvar/GRCh38/clinvar.vcf.gz
    REF=/path/to/GRCh38.fa

    loftee_dir="$VEPCACHEDIR/ensembl-vep-108.1-0/loftee-grch38/"

    invcf=$1
    outvcf=$2
    tmpvcf=${outvcf%%.vcf.gz}.tmp.vcf.gz
    threads=4

    # Run VEP annotation
    vep --assembly GRCh38 --fork $threads -i $invcf -o $tmpvcf --vcf --merged --force_overwrite --offline --use_transcript_ref --total_length --nearest symbol --compress_output bgzip \
    --everything \
    --fasta $REF \
    --dir_cache $VEPCACHEDIR \
    --dir_plugins $VEPPLUGINS \
    --plugin LoF,loftee_path:$loftee_dir,human_ancestor_fa:${loftee_dir}/GRCh38_data/human_ancestor.fa.gz,conservation_file:${loftee_dir}/GRCh38_data/loftee.sql \
    --plugin dbNSFP,$DBNSFP,ALL

    # annotate with new release clinvar
    bcftools annotate -a $CLINVAR -c CLNDISDB,CLNDN,CLNHGVS,CLNREVSTAT,CLNSIG,CLNVC,GENEINFO,MC,ORIGIN,CLNSIGCONF -Oz -o $outvcf $tmpvcf
    rm $tmpvcf
    bcftools index -f $outvcf
    
}

export -f anno

for vcf in output/pass_vcfs/*/*vcf.gz;do
    outvcf=$(sed 's/pass_vcfs/anno_vcfs/' <<< $vcf)
    mkdir -p $(dirname $outvcf)
    anno $vcf $outvcf
done
