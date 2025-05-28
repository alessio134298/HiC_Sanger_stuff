
SAMPLE_SHEET="/lustre/scratch126/gengen/teams/parts/ab77/Data/New_48336_nf_hic_samplesheet_hg38_gal6.csv"
INDEX="/lustre/scratch126/gengen/teams/parts/ab77/bwa_index/GRCh38.GRCg6a.full.renamed.merged.plusPuro.fa."
OUTPUT_FOLDER="/lustre/scratch126/gengen/teams/parts/ab77/Data/Subsample_test_HiC"
OUTPUT_FQ_SUB="${OUTPUT_FOLDER}/Fastq_subsampled"
M=128000
q="normal"
n=10

# Activate environment
module load HGI/softpack/groups/escramble/eSCRAMBLE/8

# mkdir -p ${OUTPUT_FOLDER}

# #create the 3 directories, 1 for each genome
# for p in hg38 hg38_GRCg6a GRCg6a;
#     do
#     mkdir -p ${OUTPUT_FOLDER}/${p};
# done

# ###################################################
# # Alignment, filtering, BigWig generation, idxstats 
# sed 1d ${SAMPLE_SHEET} | while IFS=',' read -r sample fastq_1 fastq_2;
#     do
#     mkdir -p ${OUTPUT_FOLDER}/hg38_GRCg6a/${sample}

#     samplename_1=$(basename "$fastq_1" | awk -F ".fastq.gz" '{print $1}')
#     samplename_2=$(basename "$fastq_2" | awk -F ".fastq.gz" '{print $1}')

#     echo "Mate 1 is: ${fastq_1}, Mate 2 is: ${fastq_2}"

#     alignment="

#     zcat ${fastq_1} | head -n 4000000 | grep @A00 | sed 's/^@//' > ${OUTPUT_FQ_SUB}/${samplename_1}.IDs.txt



    # This function is from a library that was already installed on the farm (i don't remember the github page but it was as well)
    
#     filterbyname.sh in=${fastq_1} out=${OUTPUT_FQ_SUB}/${samplename_1}.subsampled.fastq.gz names=${OUTPUT_FQ_SUB}/${samplename_1}.IDs.txt overwrite=t include=t
#     filterbyname.sh in=${fastq_2} out=${OUTPUT_FQ_SUB}/${samplename_2}.subsampled.fastq.gz names=${OUTPUT_FQ_SUB}/${samplename_1}.IDs.txt overwrite=t include=t

#     #MATE 1
#     bwa mem -t6 -A1 -B4 -E50 -L0  ${INDEX} ${OUTPUT_FQ_SUB}/${samplename_1}.subsampled.fastq.gz 2>>${OUTPUT_FOLDER}/hg38_GRCg6a/${sample}/bwa_R1.log | \
#     samtools view -hb - > ${OUTPUT_FOLDER}/hg38_GRCg6a/${sample}/${samplename_1}.sub.bam

#     samtools sort ${OUTPUT_FOLDER}/hg38_GRCg6a/${sample}/${samplename_1}.sub.bam -o \
#     ${OUTPUT_FOLDER}/hg38_GRCg6a/${sample}/${samplename_1}.sub.sorted.bam -@ 4

#     samtools index -@ 2 ${OUTPUT_FOLDER}/hg38_GRCg6a/${sample}/${samplename_1}.sub.sorted.bam

#     samtools idxstats -@ 2 ${OUTPUT_FOLDER}/hg38_GRCg6a/${sample}/${samplename_1}.sub.sorted.bam \
#     > ${OUTPUT_FOLDER}/hg38_GRCg6a/${sample}/${samplename_1}.stats

#     bamCoverage -b ${OUTPUT_FOLDER}/hg38_GRCg6a/${sample}/${samplename_1}.sub.sorted.bam -o \
#     ${OUTPUT_FOLDER}/hg38_GRCg6a/${sample}/${samplename_1}.sub.bw -p 4

#     #MATE 2
#     bwa mem -t6 -A1 -B4 -E50 -L0  ${INDEX} ${OUTPUT_FQ_SUB}/${samplename_2}.subsampled.fastq.gz 2>>${OUTPUT_FOLDER}/hg38_GRCg6a/${sample}/bwa_R2.log | \
#     samtools view -hb - > ${OUTPUT_FOLDER}/hg38_GRCg6a/${sample}/${samplename_2}.sub.bam

#     samtools sort ${OUTPUT_FOLDER}/hg38_GRCg6a/${sample}/${samplename_2}.sub.bam -o \
#     ${OUTPUT_FOLDER}/hg38_GRCg6a/${sample}/${samplename_2}.sub.sorted.bam -@ 4

#     samtools index -@ 2 ${OUTPUT_FOLDER}/hg38_GRCg6a/${sample}/${samplename_2}.sub.sorted.bam

#     samtools idxstats -@ 2 ${OUTPUT_FOLDER}/hg38_GRCg6a/${sample}/${samplename_2}.sub.sorted.bam \
#     > ${OUTPUT_FOLDER}/hg38_GRCg6a/${sample}/${samplename_2}.stats

#     bamCoverage -b ${OUTPUT_FOLDER}/hg38_GRCg6a/${sample}/${samplename_2}.sub.sorted.bam -o \
#     ${OUTPUT_FOLDER}/hg38_GRCg6a/${sample}/${samplename_2}.sub.bw -p 4
#     "

#     # Submit the job
#     bsub -n $n -J "alignment_${sample}" \
#         -o ${OUTPUT_FOLDER}/hg38_GRCg6a/${sample}/$sample.log.out \
#         -e ${OUTPUT_FOLDER}/hg38_GRCg6a/${sample}/$sample.log.err \
#         -M $M -R "select[mem>$M] rusage[mem=$M]" \
#         -q $q \
#         bash -c "$alignment"
    
#     echo "Aligning and filtering: ${sample}"

# done 
# #########################################

# #########################################
# #Splitting the files in hg38 and galGal6 sorting by name

# HG38_CHR="/lustre/scratch126/gengen/teams/parts/ab77/Data/human_chr.txt"
# GG6_CHR="/lustre/scratch126/gengen/teams/parts/ab77/Data/chicken_chr.txt"
# HG38_GG6_CHR="/lustre/scratch126/gengen/teams/parts/ab77/Data/human_chicken_chr.txt"

# sed 1d ${SAMPLE_SHEET} | while IFS=',' read -r sample fastq_1 fastq_2;
#     do
#     mkdir -p ${OUTPUT_FOLDER}/hg38/${sample}
#     mkdir -p ${OUTPUT_FOLDER}/GRCg6a/${sample}

#     samplename_1=$(basename "$fastq_1" | awk -F ".fastq.gz" '{print $1}')
#     samplename_2=$(basename "$fastq_2" | awk -F ".fastq.gz" '{print $1}')

#     echo "Mate 1 is: ${fastq_1}, Mate 2 is: ${fastq_2}"

#     #splitting
#     splitting_sorting="
#     # samtools view -hb ${OUTPUT_FOLDER}/hg38_GRCg6a/${sample}/${samplename_1}.sub.sorted.bam $(cat ${HG38_CHR}) > ${OUTPUT_FOLDER}/hg38/${sample}/${samplename_1}.hg38.sub.sorted.bam
#     # samtools view -hb ${OUTPUT_FOLDER}/hg38_GRCg6a/${sample}/${samplename_1}.sub.sorted.bam $(cat ${GG6_CHR}) > ${OUTPUT_FOLDER}/GRCg6a/${sample}/${samplename_1}.GRCg6a.sub.sorted.bam
#     # samtools view -hb ${OUTPUT_FOLDER}/hg38_GRCg6a/${sample}/${samplename_1}.sub.sorted.bam $(cat ${HG38_GG6_CHR}) > ${OUTPUT_FOLDER}/hg38_GRCg6a/${sample}/${samplename_1}.hg38_GRCg6a.sub.sorted.bam

#     # samtools view -hb ${OUTPUT_FOLDER}/hg38_GRCg6a/${sample}/${samplename_2}.sub.sorted.bam $(cat ${HG38_CHR}) > ${OUTPUT_FOLDER}/hg38/${sample}/${samplename_2}.hg38.sub.sorted.bam
#     # samtools view -hb ${OUTPUT_FOLDER}/hg38_GRCg6a/${sample}/${samplename_2}.sub.sorted.bam $(cat ${GG6_CHR}) > ${OUTPUT_FOLDER}/GRCg6a/${sample}/${samplename_2}.GRCg6a.sub.sorted.bam
#     # samtools view -hb ${OUTPUT_FOLDER}/hg38_GRCg6a/${sample}/${samplename_2}.sub.sorted.bam $(cat ${HG38_GG6_CHR}) > ${OUTPUT_FOLDER}/hg38_GRCg6a/${sample}/${samplename_2}.hg38_GRCg6a.sub.sorted.bam

#     # #sort by name
#     # samtools sort -n ${OUTPUT_FOLDER}/hg38/${sample}/${samplename_1}.hg38.sub.sorted.bam -o ${OUTPUT_FOLDER}/hg38/${sample}/${samplename_1}.hg38.sub.sortedbyname.bam
#     # samtools sort -n ${OUTPUT_FOLDER}/GRCg6a/${sample}/${samplename_1}.GRCg6a.sub.sorted.bam -o ${OUTPUT_FOLDER}/GRCg6a/${sample}/${samplename_1}.GRCg6a.sub.sortedbyname.bam
#     # samtools sort -n ${OUTPUT_FOLDER}/hg38_GRCg6a/${sample}/${samplename_1}.hg38_GRCg6a.sub.sorted.bam -o ${OUTPUT_FOLDER}/hg38_GRCg6a/${sample}/${samplename_1}.hg38_GRCg6a.sub.sortedbyname.bam

#     # samtools sort -n ${OUTPUT_FOLDER}/hg38/${sample}/${samplename_2}.hg38.sub.sorted.bam -o ${OUTPUT_FOLDER}/hg38/${sample}/${samplename_2}.hg38.sub.sortedbyname.bam
#     # samtools sort -n ${OUTPUT_FOLDER}/GRCg6a/${sample}/${samplename_2}.GRCg6a.sub.sorted.bam -o ${OUTPUT_FOLDER}/GRCg6a/${sample}/${samplename_2}.GRCg6a.sub.sortedbyname.bam
#     # samtools sort -n ${OUTPUT_FOLDER}/hg38_GRCg6a/${sample}/${samplename_2}.hg38_GRCg6a.sub.sorted.bam -o ${OUTPUT_FOLDER}/hg38_GRCg6a/${sample}/${samplename_2}.hg38_GRCg6a.sub.sortedbyname.bam

#     samtools sort -n ${OUTPUT_FOLDER}/hg38_GRCg6a/${sample}/${samplename_1}.sub.bam -o ${OUTPUT_FOLDER}/hg38_GRCg6a/${sample}/${samplename_1}.hg38_GRCg6a.sub.sortedbyname.bam
#     samtools sort -n ${OUTPUT_FOLDER}/hg38_GRCg6a/${sample}/${samplename_2}.sub.bam -o ${OUTPUT_FOLDER}/hg38_GRCg6a/${sample}/${samplename_2}.hg38_GRCg6a.sub.sortedbyname.bam
#     "

#     # Submit the job
#     bsub -n $n -J "filtering_${sample}" \
#         -w "done(alignment_${sample})" \
#         -o ${OUTPUT_FOLDER}/hg38_GRCg6a/${sample}/$sample.log.out \
#         -e ${OUTPUT_FOLDER}/hg38_GRCg6a/${sample}/$sample.log.err \
#         -M $M -R "select[mem>$M] rusage[mem=$M]" \
#         -q $q \
#         bash -c "$splitting_sorting";
# done

#######################################

#######################################
#HiCExplorer

REST_SITES_DIR="/lustre/scratch126/gengen/teams/parts/ab77/Data/RestSites"

# mkdir -p ${REST_SITES_DIR}

# findrestsite="
# # DpnII (4 bp) - GATC
# hicFindRestSite \
# 	--fasta /lustre/scratch126/gengen/projects/graft/Dataset/reference/hg38_galGal6_full_plusPuro/fasta/GRCh38.GRCg6a.full.renamed.merged.plusPuro.fa \
# 	--searchPattern GATC \
# 	-o  ${REST_SITES_DIR}/DpnII_site_positions.bed 
# # HinfI (5 bp) - GANTC
# hicFindRestSite \
# 	--fasta /lustre/scratch126/gengen/projects/graft/Dataset/reference/hg38_galGal6_full_plusPuro/fasta/GRCh38.GRCg6a.full.renamed.merged.plusPuro.fa \
# 	--searchPattern GA.TC \
# 	-o  ${REST_SITES_DIR}/HinfI_site_positions.bed
# # Combine rest sites
# cat  ${REST_SITES_DIR}/DpnII_site_positions.bed  ${REST_SITES_DIR}/HinfI_site_positions.bed >  ${REST_SITES_DIR}/Arima_site_positions.bed
# sort -k1,1 -k2,2n  ${REST_SITES_DIR}/Arima_site_positions.bed > ${REST_SITES_DIR}/Arima_site_positions.sorted.bed
# "

# bsub -n 1 -J "findrestsite" \
#     -w "done()"\
#     -o ${REST_SITES_DIR}/findrestsite.log.out \
#     -e ${REST_SITES_DIR}/findrestsite.log.err \
#     -M $M -R "select[mem>$M] rusage[mem=$M]" \
#     -q $q \
#     bash -c "${findrestsite}"

# Build contact matrix
genome="hg38_GRCg6a"

sed 1d ${SAMPLE_SHEET} | while IFS="," read -r sample fastq_1 fastq_2;
    do
    samplename_1=$(basename "$fastq_1" | awk -F ".fastq.gz" '{print $1}')
    samplename_2=$(basename "$fastq_2" | awk -F ".fastq.gz" '{print $1}')
    
    build_matrix="
    # hicBuildMatrix --samFiles ${OUTPUT_FOLDER}/${genome}/${sample}/${samplename_1}.${genome}.sub.sortedbyname.bam \
    #             ${OUTPUT_FOLDER}/${genome}/${sample}/${samplename_2}.${genome}.sub.sortedbyname.bam \
    #             --outBam ${OUTPUT_FOLDER}/${genome}/${sample}/${sample}.${genome}.HiC_1kb.bam \
    #             --binSize 1000 \
    #             --restrictionSequence GATC GA.TC \
    #             --danglingSequence GATC A.T \
    #             --restrictionCutFile ${REST_SITES_DIR}/Arima_site_positions.sorted.bed \
    #             --threads 6 \
    #             --inputBufferSize 400000 \
    #             -o ${OUTPUT_FOLDER}/${genome}/${sample}/${sample}.${genome}.HiC_1kb.h5 \
    #             --QCfolder ${OUTPUT_FOLDER}/${genome}/${sample}/${sample}.${genome}.HiC_1kb_QC
    
    samtools sort ${OUTPUT_FOLDER}/${genome}/${sample}/${sample}.${genome}.HiC_1kb.bam -o \
    ${OUTPUT_FOLDER}/${genome}/${sample}/${sample}.${genome}.HiC_1kb.sorted.bam
    samtools index ${OUTPUT_FOLDER}/${genome}/${sample}/${sample}.${genome}.HiC_1kb.sorted.bam
    "
    # submit job only when the previous has finished
    bsub -n $n -J "hicBuildMatrix_${sample}" \
        -o ${OUTPUT_FOLDER}/${genome}/${sample}/${sample}.log.out \
        -e ${OUTPUT_FOLDER}/${genome}/${sample}/${sample}.log.err \
        -M $M -R "select[mem>$M] rusage[mem=$M]" \
        -q $q \
        bash -c "${build_matrix}"
done