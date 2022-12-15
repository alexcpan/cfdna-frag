import argparse
import pandas as pd
import csv
import sys

parser = argparse.ArgumentParser()
parser.add_argument('--bed_file',help='input bedfile',required=True)
parser.add_argument('--distro_file',help='output distribution file',required=True)
args = parser.parse_args()

bed_file = args.bed_file
distro_file = args.distro_file

print("bed_file",bed_file)
print("distro_file",distro_file)

# df = (pd.read_csv(bed_file,sep="\t",usecols=[3], header=None, names=["gc"])
#     .round(2)
#     .value_counts().rename_axis('gc_strata').reset_index(name='n')
# )

# df["fract_frags"] = df["n"]/df["n"].sum()
# df.insert(len(df.columns),"fract_frags",df["n"]/df["n"].sum())
# lib_id = bed_file.split("/")[-1].split("_")[0]
# df.insert(0,"library_id",lib_id)
# df.sort_values(by=['gc_strata'],inplace=True)
# df.to_csv(distro_file, index=False,columns=["library_id","gc_strata","fract_frags"])


#### NOTES:
#
# Expect the 5 columns in bedfile to be:
# 0_chr 1_start 2_end 3_pct_gc 4_seq_len
# we only care about pct_gc
#
# We have to read in chunks otherwise we get memory issues with the massive bams.
# No worries about saving results in df because the size is bounded (100 gc_strata).
#

df = pd.DataFrame(columns = ["gc_strata","n"])
CHUNKSIZE = 10**6 #(1mb chunks)
i = 1
for chunk in pd.read_csv(bed_file,sep="\t",usecols=[3], header=None, names=["gc"], chunksize=CHUNKSIZE):
    if i % 100 == 0:
        print("chunk number: ", i,flush=True)
    chunk_df = chunk.round(2).value_counts().rename_axis('gc_strata').reset_index(name='n')
    df = pd.concat([df, chunk_df]).groupby(["gc_strata"]).sum().reset_index()
    i += 1

(df.assign(fract_frags=df["n"]/df["n"].sum())
    .assign(library_id = bed_file.split("/")[-1].split("_")[0])
    .sort_values(by=['gc_strata'])
    .to_csv(distro_file, index=False,columns=["library_id","gc_strata","fract_frags"],quoting=csv.QUOTE_NONNUMERIC)
)
