#!/bin/bash
SAMPLE_SHEET="/lustre/scratch126/gengen/projects/graft/Analysis/ab77/hic_analysis/new_data/galGal6_human_samples.csv"
FASTA_HG38_GALGAL6="/lustre/scratch126/gengen/projects/graft/Dataset/reference/hg38_galGal6_full_plusPuro/fasta/GRCh38.GRCg6a.full.renamed.merged.plusPuro.fa"
FASTA_HG38_MM39="/lustre/scratch126/gengen/projects/graft/Dataset/reference/hg38_mm39/fasta/hg38_mm39_renamed.fa"
INDEX_HG38_GALGAL6="/lustre/scratch126/gengen/teams/parts/ab77/bwa_index/GRCh38.GRCg6a.full.renamed.merged.plusPuro.fa."
INDEX_HG38_MM39="/lustre/scratch126/gengen/projects/graft/Dataset/reference/hg38_mm39/bwa_index/hg38_mm39_renamed.fa"
FOLDER="/lustre/scratch126/gengen/projects/graft/Analysis/ab77/hic_analysis/new_data/Processed_samples/chicken"
REST_SITES_DIR="/lustre/scratch126/gengen/projects/graft/Analysis/ab77/hic_analysis/RestSites"
CHROMSIZES_HG38_GALGAL6="/lustre/scratch126/gengen/projects/graft/Dataset/reference/hg38_galGal6_full_plusPuro/GRCh38.GRCg6a.full.renamed.merged.plusPuro.chromsizes.txt"

module load HGI/softpack/groups/escramble/eSCRAMBLE/8

M=128000
q="long"
n=20

while IFS=',' read -r sample fastq_1 fastq_2 genome
    do

    mkdir -p ${FOLDER}/${sample}/Akita
    mkdir -p ${FOLDER}/${sample}/Akita/logs

    samplename_1=$(basename "$fastq_1" | awk -F ".fastq.gz" '{print $1}')
    samplename_2=$(basename "$fastq_2" | awk -F ".fastq.gz" '{print $1}')

    echo "Processing: ${sample}, Mate 1 is: ${fastq_1}, Mate 2 is: ${fastq_2}, genome is: ${genome}"

    if [ "${genome}" == "hg38_galGal6" ]; then
        INDEX=${INDEX_HG38_GALGAL6}
    elif [ "${genome}" == "hg38_mm39" ]; then
        INDEX=${INDEX_HG38_MM39}
    else
        echo "error: $genome not recognized"
        exit 1
    fi

    Akita_preprocess_1="
    # Alignment, sortingbyname

    bwa mem -t20 -A1 -B4 -E50 -L0  ${INDEX} ${fastq_1} 2>>${FOLDER}/${sample}/Akita/logs/bwa_R1.log | \
    samtools view -hb - | \
    samtools sort -n - -o ${FOLDER}/${sample}/Akita/${samplename_1}.sortedbyname.bam \
    -@ 20 2>>${FOLDER}/${sample}/Akita/logs/sortbyname_R1.log

    bwa mem -t20 -A1 -B4 -E50 -L0  ${INDEX} ${fastq_2} 2>>${FOLDER}/${sample}/Akita/logs/bwa_R2.log | \
    samtools view -hb - | \
    samtools sort -n - -o ${FOLDER}/${sample}/Akita/${samplename_2}.sortedbyname.bam \
    -@ 20 2>>${FOLDER}/${sample}/Akita/logs/sortbyname_R2.log
    "

    bsub -n $n -J "Alignment_sort_${sample}" \
        -o ${FOLDER}/${sample}/Akita/logs/Align.log.out \
        -e ${FOLDER}/${sample}/Akita/logs/Align.log.err \
        -M $M -R "select[mem>$M] rusage[mem=$M]" \
        -q $q \
        bash -c "${Akita_preprocess_1}"

done < "$SAMPLE_SHEET"


M=256000
q="long"
n=8


while IFS=',' read -r sample fastq_1 fastq_2 genome;
    do
    
    samplename_1=$(basename "$fastq_1" | awk -F ".fastq.gz" '{print $1}')
    samplename_2=$(basename "$fastq_2" | awk -F ".fastq.gz" '{print $1}')

    # HiCExplorer

    Res=2048

    HiC_Matrix_1024="
    hicBuildMatrix --samFiles ${FOLDER}/${sample}/Akita/${samplename_1}.sortedbyname.bam \
                ${FOLDER}/${sample}/Akita/${samplename_2}.sortedbyname.bam \
                --binSize 1024 \
                --restrictionSequence GATC GA.TC \
                --danglingSequence GATC A.T \
                --restrictionCutFile ${REST_SITES_DIR}/Arima_site_positions_${genome}.sorted.bed \
                --threads 12 \
                --inputBufferSize 400000 \
                -o ${FOLDER}/${sample}/Akita/${sample}.${genome}.HiC_1024b.h5 \
                2>>${FOLDER}/${sample}/Akita/logs/buildmatrix.log

    # rm ${FOLDER}/${sample}/Akita/${samplename_1}.sortedbyname.bam ${FOLDER}/${sample}/Akita/${samplename_2}.sortedbyname.bam

    hicConvertFormat \
		-m ${FOLDER}/${sample}/Akita/${sample}.${genome}.HiC_1024b.h5 \
		--inputFormat h5 \
		--outputFormat mcool \
		-r 1024 2048 5120 10240 20480 51200 102400 204800 512000 1024000 \
		-o ${FOLDER}/${sample}/Akita/${sample}.${genome}.HiC_1024b.mcool
    
    
    cooler cp ${FOLDER}/${sample}/Akita/${sample}.${genome}.HiC_1024b.mcool::/resolutions/${Res} \
    ${FOLDER}/${sample}/Akita/${sample}.${genome}.HiC_${Res}b.cool

    cp ${FOLDER}/${sample}/Akita/${sample}.${genome}.HiC_${Res}b.cool \
    ${FOLDER}/${sample}/Akita/${sample}.${genome}.HiC_${Res}b.balanced.cool

    cooler balance -p 8 ${FOLDER}/${sample}/Akita/${sample}.${genome}.HiC_${Res}b.balanced.cool --force
    "

    if bjobs -w | grep -q "Alignment_sort_${sample}"; then
        wait_condition="-w done(Alignment_sort_${sample})"
    else
        wait_condition=""
    fi

    # Submit the job
    bsub -n $n -J "HiC_Matrix_1024_${sample}" \
        $wait_condition \
        -o ${FOLDER}/${sample}/Akita/${sample}_Matrix.log.out \
        -e ${FOLDER}/${sample}/Akita/${sample}_Matrix.log.err \
        -M $M -R "select[mem>$M] rusage[mem=$M]" \
        -q $q \
        bash -c "${HiC_Matrix_1024}"

done < ${SAMPLE_SHEET}