#!/bin/bash
SAMPLE_SHEET="/lustre/scratch126/gengen/projects/graft/Analysis/ab77/hic_analysis/new_data/galGal6_human_samples.csv"
OUTPUT_FOLDER="/lustre/scratch126/gengen/projects/graft/Analysis/ab77/hic_analysis/new_data/Processed_samples/chicken"
M=128000
q="normal"
n=20

source /software/hgi/installs/conda-audited/miniforge/etc/profile.d/conda.sh
conda activate hic_new

Res=500000

while IFS=',' read -r sample fastq_1 fastq_2 genome;
    do
    
    mkdir -p ${OUTPUT_FOLDER}/${sample}/cooler

    balance_dump="

    cooler cp ${OUTPUT_FOLDER}/${sample}/${sample}.${genome}.HiC_1kb.mcool::/resolutions/${Res} \
    ${OUTPUT_FOLDER}/${sample}/cooler/${sample}.${genome}.${Res}.cool

    cp ${OUTPUT_FOLDER}/${sample}/cooler/${sample}.${genome}.${Res}.cool \
    ${OUTPUT_FOLDER}/${sample}/cooler/${sample}.${genome}.${Res}.balanced.cool

    cp ${OUTPUT_FOLDER}/${sample}/cooler/${sample}.${genome}.${Res}.cool \
    ${OUTPUT_FOLDER}/${sample}/cooler/${sample}.${genome}.${Res}.balanced_intra.cool

    cp ${OUTPUT_FOLDER}/${sample}/cooler/${sample}.${genome}.${Res}.cool \
    ${OUTPUT_FOLDER}/${sample}/cooler/${sample}.${genome}.${Res}.balanced_inter.cool

    Ho dimenticato qua di aggiungere --cis-only e --trans-only

    cooler balance -p 30 ${OUTPUT_FOLDER}/${sample}/cooler/${sample}.${genome}.${Res}.balanced.cool --force
    cooler balance -p 30 ${OUTPUT_FOLDER}/${sample}/cooler/${sample}.${genome}.${Res}.balanced_intra.cool --force --cis-only
    cooler balance -p 30 ${OUTPUT_FOLDER}/${sample}/cooler/${sample}.${genome}.${Res}.balanced_inter.cool --force --trans-only

    cooler dump --balanced ${OUTPUT_FOLDER}/${sample}/cooler/${sample}.${genome}.${Res}.balanced.cool -o \
    ${OUTPUT_FOLDER}/${sample}/cooler/${sample}.${genome}.${Res}.balanced.txt -t pixels --join --header

    cooler dump --balanced ${OUTPUT_FOLDER}/${sample}/cooler/${sample}.${genome}.${Res}.balanced_intra.cool -o \
    ${OUTPUT_FOLDER}/${sample}/cooler/${sample}.${genome}.${Res}.balanced_intra.txt -t pixels --join --header

    cooler dump --balanced ${OUTPUT_FOLDER}/${sample}/cooler/${sample}.${genome}.${Res}.balanced_inter.cool -o \
    ${OUTPUT_FOLDER}/${sample}/cooler/${sample}.${genome}.${Res}.balanced_inter.txt -t pixels --join --header
    "

    bsub -n $n -J "Cooler_balance_dump_${sample}" \
        -o ${OUTPUT_FOLDER}/${sample}/cooler/${sample}.log.out \
        -e ${OUTPUT_FOLDER}/${sample}/cooler/${sample}.log.err \
        -M $M -R "select[mem>$M] rusage[mem=$M]" \
        -q $q \
        bash -c "${balance_dump}"

done < ${SAMPLE_SHEET}