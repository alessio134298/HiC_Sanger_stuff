import cooler # type: ignore
import numpy as np

Matrix = "/lustre/scratch126/gengen/projects/graft/Analysis/ab77/hic_analysis/new_data/Processed_samples/5c3.rep1/5c3.rep1.hg38_galGal6.HiC_1kb.corrected.mcool"

# with 100000 worked trying with 1000
c = cooler.Cooler(f"{Matrix}::resolutions/1000")

# Balance does not work, probably because the matrix was not originary made with cooler, so it has not been "balanced" with the tool
mat = c.matrix(balance=False)[:]
np.save("/lustre/scratch126/gengen/projects/graft/Analysis/ab77/hic_analysis/new_data/Processed_samples/5c3.rep1/5c3.rep1.hg38_galGal6.HiC_1kb.corrected_1000.npy", mat)


# bsub -n 1 -J "npy" -e npy.log.err -o npy.log.out  -M 256000 -R "select[mem>256000] rusage[mem=256000]" -q "normal" bash -c "python3 /lustre/scratch126/gengen/teams/parts/ab77/scripts/trans-C/convert_to_numpy.py"