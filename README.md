# CHN_ECS

This repository contains the code used for [Populational pan-ethnic screening panel enabled by deep whole genome sequencing]().

## Obain public datasets

Three public datasets were involed in the study:

- The China Metabolic Analytics Project (ChinaMAP)

    ChinaMap provides a site-only VCF for downloading, after sign in [mbiobank](http://www.mbiobank.com/), the VCF file can be downloaded from [this page](http://www.mbiobank.com/download/)

- Westlake BioBank for Chinese (WBBC)

    WBBC also provides sites and allele frequencies VCFs. As we were focus on AR disease, only autosomal VCFs of GRCh38 were downloaded from [this page](https://wbbc.westlake.edu.cn/downloads.html)

- The Genome Aggregation Database (gnomAD)

    The variant dataset files of the latest release of gnomAD (v3) were avaliable from [this page](http://www.gnomad-sg.org/downloads#v3-variants). They contain information of all subsets: non-cancer, non-neuro, non-v2, non-TOPMed, controls/biobanks, 1KG, and HGDP. The non-cancer subset was later extracted for analysis in this study.

## Content description

### 0. **`data/`**

- `candidate_gene_excludeXY.list`: a list of candidate AR genes
- `1-s2.0-S1534580722004506-mmc2.xls`: supplementary data of [a paper](https://doi.org/10.1016/j.devcel.2022.06.010)

### 1. **`scripts/`**

- `parse_vcf_by_gene.py`: This script is used to parse variants by given a list of gene SYMBOL. It takes a VCF annotated by VEP/clinvar and a gene list files as input, then extract variants by looking at the SYMBOL or GENEINFO of the INFO fields if they appears in the given gene list.

- `screen_pathogenic_variant.py`: This script is used to identify pathogenic variants according to the method describe in the paper. a estimated GCR is also calculate of each gene.

- `llps_annotation.py`: This script is used to check if a pathogenic variants is reported in [Banani SF, Afeyan LK, Hawken SW, Henninger JE, Dall'Agnese A, Clark VE, Platt JM, Oksuz O, Hannett NM, Sagi I, Lee TI, Young RA. Genetic variation associated with condensate dysregulation in disease. Dev Cell. 2022 Jul 25;57(14):1776-1788.e8.](https://doi.org/10.1016/j.devcel.2022.06.010)

### 2. **`analysis/`**

- `s0.preprocessing.sh`: Step0: Pre-processing and QC of downloaded VCFs
- `s1.variant_annotation.sh`: Step1: Annotate variants by VEP and Clinvar Database
- `s2.extract_vars_on_candidate_gene.sh`: Step2: Extract variants of the candidate genes
- `s3.screen_pathogenetic_variants_estimate_GCR.sh`: Step3: identify pathogenic variants and calculate GCR
- `s4.llps_analysis.sh`: Step4: LLPS related variants/Gene analaysis

## Requirements

* Python >= 3.7 (cyvcf2, pandas)
* VEP 108
* bcftools >= 1.16
