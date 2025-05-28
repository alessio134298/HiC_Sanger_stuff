#!/bin/bash
SAMPLE_SHEET="/lustre/scratch126/gengen/projects/graft/Analysis/ab77/hic_analysis/new_data/Samplesheet_mouse_matrix.csv"
FOLDER="/lustre/scratch126/gengen/projects/graft/Analysis/ab77/hic_analysis/new_data/Merged_Matrix_mouse"
CHROMSIZES_HG38_MM39="/lustre/scratch126/gengen/projects/graft/Dataset/reference/hg38_mm39/hg38_mm39_renamed_chromsizes.txt"
M=256000
q="long"
n=12

module load HGI/softpack/groups/escramble/eSCRAMBLE/8


CHR_mESC_Hchr21=$(echo mm39_{1..19} mm39_X mm39_Y hg38_21)
CHR_mESC_Hchr4=$(echo mm39_{1..19} mm39_X mm39_Y hg38_4)
CHR_HUMAN=$(echo hg38_{1..22} hg38_X hg38_Y)

while IFS=',' read -r sample Rep1 Rep2 Rep3 Rep4 genome

    do

    if [ "${sample}" == "WT.RPE1" ]; then
        CHR=${CHR_HUMAN}
        LOW_LIM=-1.8
        MAX_LIM=3
    elif [ "${sample}" == "RPE1.PB_NEO" ]; then
        CHR=${CHR_HUMAN}
        LOW_LIM=-2.2
        MAX_LIM=3
    elif [ "${sample}" == "CT2_1" ]; then
        CHR=${CHR_mESC_Hchr21}
        LOW_LIM=-1.6
        MAX_LIM=2.6
    elif [ "${sample}" == "CL3" ]; then
        CHR=${CHR_mESC_Hchr4}
        LOW_LIM=-1.6
        MAX_LIM=2
    elif [ "${sample}" == "7M5_1_2" ]; then
        CHR=${CHR_HUMAN}
        LOW_LIM=-1.8
        MAX_LIM=2
    elif [ "${sample}" == "2A8" ]; then
        CHR=${CHR_HUMAN}
        LOW_LIM=-1.2
        MAX_LIM=2
    elif [ "${sample}" == "el_Cl35" ]; then
        CHR=${CHR_HUMAN}
        LOW_LIM=-2.4
        MAX_LIM=3
    elif [ "${sample}" == "el_Cl2" ]; then
        CHR=${CHR_HUMAN}
        LOW_LIM=-2.6
        MAX_LIM=3
    else 
        echo "not in the list of chromosomes"
        exit 1
    fi

    param=(
    "${LOW_LIM} ${MAX_LIM}"
    )

    # Correct ICE


    name=$(echo ${param} | tr " " "_") 

    correct_ICE="
    hicCorrectMatrix correct -m ${FOLDER}/${sample}/${sample}.${genome}.HiC_1kb.h5 \
    -o ${FOLDER}/${sample}/${sample}.${genome}.HiC_1kb.corrected.ICE_${name}.h5 \
    --correctionMethod ICE \
    --chromosomes ${CHR} \
    --filterThreshold ${param}

    bash /lustre/scratch126/gengen/teams/parts/ab77/scripts/Convert_to_hic.sh \
    ${FOLDER}/${sample}/${sample}.${genome}.HiC_1kb.corrected.ICE_${name}.h5 \
    ${CHROMSIZES_HG38_MM39} \
    h5
    "

    # echo ${correct_ICE}

    bsub -J "correct_ICE_${sample}_${name}" \
    -o ${FOLDER}/${sample}/${sample}.correct_Mat_ICE_${name}.log.out \
    -e ${FOLDER}/${sample}/${sample}.correct_Mat_ICE_${name}.log.err \
    -M $M -R "select[mem>$M] rusage[mem=$M]" \
    -q $q \
    bash -c "${correct_ICE}"

done < ${SAMPLE_SHEET}