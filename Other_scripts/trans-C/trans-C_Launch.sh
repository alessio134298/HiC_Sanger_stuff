#!bin/bash
M=128000
n=12
q="normal"

source /software/hgi/installs/conda-audited/miniforge/etc/profile.d/conda.sh
conda activate trans-C

MATRIX="/lustre/scratch126/gengen/projects/graft/Analysis/ab77/hic_analysis/new_data/Processed_samples/5c3.rep1/5c3.rep1.hg38_galGal6.HiC_1kb.corrected.npy"
CHROMSIZES="/lustre/scratch126/gengen/projects/graft/Dataset/reference/hg38_galGal6_full_plusPuro/GRCh38.GRCg6a.full.renamed.merged.plusPuro.chromsizes.txt"
BIN=100000
SEED_FILE="/lustre/scratch126/gengen/teams/parts/ab77/trans-C/gg6_chr1_400bp_windows_20subset.bed"
OUTDIR="/lustre/scratch126/gengen/teams/parts/ab77/trans-C/5c3.rep1_test2"

module load HGI/softpack/groups/escramble/eSCRAMBLE/8

mkdir -p ${OUTDIR}

cmd="
python3 /lustre/scratch126/gengen/teams/parts/ab77/trans-C/code/trans_C.py ${MATRIX} ${CHROMSIZES} ${BIN} ${SEED_FILE} ${OUTDIR}
"

bsub -n $n -J "trans-C" \
        -o ${OUTDIR}/test.log.out \
        -e ${OUTDIR}/test.log.err \
        -M $M -R "select[mem>$M] rusage[mem=$M]" \
        -q $q \
        bash -c "${cmd}"