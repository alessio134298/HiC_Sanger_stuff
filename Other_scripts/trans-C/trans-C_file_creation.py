# this script is to create a bed file to use as 'seed_file.bed" for trans-C analyis
# the tool require a bed file to use a seed regions which dtermine trans inter-chromosomal contacts
# since we want to determine the position of gg6_1 this could be an approcah to make this chromosome by itself in spatial relation  with the others

import os
import pandas as pd

# count the length of the chromosome
file = "/lustre/scratch126/gengen/projects/graft/Dataset/reference/hg38_galGal6_full_plusPuro/GRCh38.GRCg6a.full.renamed.merged.plusPuro.chromsizes.txt"

chromsizes = pd.read_csv(file, sep = '\t', header = None)
gg6_1 = chromsizes.loc[chromsizes[0] == "gg6_1"]
length = gg6_1.iloc[0,1].item()
print(f"Length of the chrmsome is: {length}")
# Length of the chrmsome is: 197608386

# create the dataframe
binned_bed = pd.DataFrame(columns = ["chr", "start", "end"])

dict_list = []

for i in range(1, length, 400):
    row_dict = {'chr': "gg6_1", 'start':i, 'end':i  + 399}
    dict_list.append(row_dict)

binned_bed = pd.DataFrame.from_dict(dict_list)

binned_bed.drop(binned_bed.tail(1).index,inplace=True)

binned_bed.to_csv("/lustre/scratch126/gengen/teams/parts/ab77/trans-C/gg6_chr1_400bp_windows.bed", sep='\t', index=False, header=False)
