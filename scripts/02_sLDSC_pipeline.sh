#!/bin/bash

eval "$(~/anaconda3/bin/conda shell.bash hook)"
source activate ldsc

#paths to the directories
LDSC="/path/to/ldsc_software"
REQUIRED_FILES="/path/to/ldsc_required_files"
ANNOT_FILES="/path/to/annotation_files"
HERITABILITY_OUTPUT="/path/to/heritability_output"
BEDFILES="/path/to/bedfiles"
DISEASE="/path/to/gwas_sumstats"

# stage-specific consensus bed files, one annotation per stage
stage=()
for file in "$BEDFILES/"*_consensus.bed; do
    if [ -e "$file" ]; then
        filename=$(basename "$file")
        name="${filename%_consensus.bed}"
        if [[ ! "$name" =~ _ac_hg19$ ]]; then
            stage+=("$name")
        fi
    fi
done

chrom=({1..22})

## Step 1: Creating an annot file
parallel --env LDSC --env BEDFILES --env REQUIRED_FILES --env ANNOT_FILES --env HERITABILITY_OUTPUT -j 20 "echo \"Running ldsc for: {1} and {2} \"; \
python ${LDSC}/make_annot.py \
  --bed-file ${BEDFILES}/{1}_consensus.bed \
  --bimfile ${REQUIRED_FILES}/1000G_EUR_Phase3_plink/1000G.EUR.QC.{2}.bim \
  --annot-file ${ANNOT_FILES}/olig2_ageing.{1}.{2}.annot.gz" ::: "${stage[@]}" ::: "${chrom[@]}"

## Step 2: Generating ldsc scores
parallel --env LDSC --env BEDFILES --env REQUIRED_FILES --env ANNOT_FILES --env HERITABILITY_OUTPUT -j 20 "echo \"Running ldsc for: {1} and {2} \"; \
python ${LDSC}/ldsc.py \
  --l2 \
  --bfile ${REQUIRED_FILES}/1000G_EUR_Phase3_plink/1000G.EUR.QC.{2} \
  --ld-wind-cm 1 \
  --annot ${ANNOT_FILES}/olig2_ageing.{1}.{2}.annot.gz  \
  --thin-annot \
  --out ${ANNOT_FILES}/olig2_ageing.{1}.{2} \
  --print-snps ${REQUIRED_FILES}/list.txt" ::: "${stage[@]}" ::: "${chrom[@]}"

## Step 3: Partitioned heritability regression
export diseases=("AD_Jansen2019" "AD_Kunkle2019" "AD_Bellenguez2022" "SVD_Sargurupremraj2022" "PVS_Duperron2023" "PD_GP22025" "MS_Andlauer2016" "ALS_Rheenen2021" "SCZ_Trubetskoy2022" "Covid_Castineira2023" "Longevity_Deelan2019_90th")
type=("fetal" "pediatric" "adult" "aged")

export -f python

parallel --env LDSC --env DISEASE --env REQUIRED_FILES --env ANNOT_FILES --env HERITABILITY_OUTPUT -j 20 \
"echo '>>> Running LDSC for {1} with {2}'; \
python ${LDSC}/ldsc.py \
  --h2 ${DISEASE}/{1}_ldsc.sumstats.gz \
  --ref-ld-chr ${REQUIRED_FILES}/1000G_Phase3_baselineLD_v2.2_ldscores/baselineLD.,${ANNOT_FILES}/olig2_ageing.{2}. \
  --out ${HERITABILITY_OUTPUT}/olig2_ageing.{2}.{1} \
  --overlap-annot \
  --frqfile-chr ${REQUIRED_FILES}/1000G_Phase3_frq/1000G.EUR.QC. \
  --w-ld-chr ${REQUIRED_FILES}/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
  --print-coefficients" \
::: "${diseases[@]}" ::: "${type[@]}"
