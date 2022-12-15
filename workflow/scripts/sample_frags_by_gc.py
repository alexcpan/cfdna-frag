import argparse
import pandas as pd
import csv
import sys
import os
import numpy as np
import time
import shutil
import glob

parser = argparse.ArgumentParser()
parser.add_argument('--healthy_csv',help='healthy gc strata csv',required=True)
parser.add_argument('--bed_file',help='output distribution file',required=True)
parser.add_argument('--sampled_bed',help='output sampled bed file',required=True)
args = parser.parse_args()

healthy_csv = args.healthy_csv
bed_file    = args.bed_file
sampled_bed = args.sampled_bed

# TODO: add checks for healthy_df/bedfile/output formatting

# Expect 2 columns in healthy df:
# 0_gc_strata,1_med_frag_fract

# Expect the 5 columns in bedfile to be:
# 0_chr 1_start 2_end 3_pct_gc 4_seq_len
# we only care about pct_gc

print('\n------------------ Params and globals------------------')
print("healthy_csv:",healthy_csv)
print("bed_file   :",bed_file)
print("sampled_bed:",sampled_bed)

LIB_ID = bed_file.split("/")[-1].split("_")[0]
TEMP_DIR = bed_file.split(LIB_ID)[0] + LIB_ID + "_tmp"
print("lib_id     :",LIB_ID)
print("temp_dir   :",TEMP_DIR)
if not os.path.exists(TEMP_DIR):
    os.makedirs(TEMP_DIR)

RS = 0
np.random.seed(RS)
print("rand_seed  :",TEMP_DIR)

print("\n------------------------------------")

###############################################################################
#### TESTING
"""
# test
cd ~/data/cfdna-frag/workflow/scripts/
python3 sample_frags_by_gc.py \
--healthy_csv /data/panalc/cfdna-frag/test/distro/healthy_med.csv \
--bed_file /data/panalc/cfdna-frag/test/frag/lib001_frag.bed \
--sampled_bed /data/panalc/cfdna-frag/test/frag/lib001_frag_sampled.2.bed

# lib250
python3 sample_frags_by_gc.py \
--healthy_csv /data/panalc/cfdna-frag/MPNST_frag/distro/healthy_med.csv \
--bed_file /data/panalc/cfdna-frag/MPNST_frag/frag/lib250_frag.bed \
--sampled_bed /data/panalc/cfdna-frag/MPNST_frag/frag/lib250_frag_sampled.bed
"""
###############################################################################
#### NOTES:
# This basically ends up doing two levels of stratified sampling.
# More information:
# https://stats.stackexchange.com/questions/185229/when-should-you-choose-stratified-sampling-over-random-sampling

# Overall, this is IO-heavy and kind of slow (takes a few minutes)
# but it provides a cap on memory usage. Some speedup via
# multi-threading is possible and was tested, but without a priori
# knowledge of the memory limits, it's still possible to get OOM-ed
# since we won't know how many processes (and thus how much memory)
# we'll be able to allocate.

# In theory we should be able to get some of this info from the "resource"
# package: https://docs.python.org/3/library/resource.html
# But this already runs in reasonable time (on the order of minutes not hours)
# so it might not be worth optimizing further at this point.

###############################################################################
#### SETUP DIST DF
# load healthy gc strata and set gc strata (2 decimal) as index
# !!!: need to round bc of floating point nonsense
STRATA = [round(s,2) for s in np.arange(0.0, 1.01, 0.01)]
healthy_df = pd.DataFrame(index=STRATA)
healthy_df = healthy_df.join(pd.read_csv(healthy_csv, index_col=0)).fillna(0)
# healthy_df = pd.read_csv(healthy_csv, index_col=0)


###############################################################################
#### SPLIT BED INTO STRATA FILES
COLS =["chr","start","end","gc","seq_len"]
CHUNKSIZE = 10**6 #(1mb chunks)

t0 = time.time()
print('------------------ Splitting .bed into strata ------------------\n')

chunk_i = 0
total_rows = 0
for chunk in pd.read_csv(bed_file,
                        sep="\t", 
                        header=None, 
                        names=COLS, 
                        chunksize=CHUNKSIZE):
    # if chunk_i >= MAX_CHUNKS:
    #     break
    print("chunk number: ", chunk_i,flush=True)
    chunk_i += 1
    chunk_df = chunk.round({"gc":2})
    chunk_df["gc_weight"] = chunk_df["gc"].map(healthy_df["med_frag_fract"])
    rows = 0
    for strata in np.arange(0.0, 1.01, 0.01):
        strata = round(strata,2) # (!!! ensure floats are properly rounded. will drop rows if removed)
        strata_df = chunk_df.loc[chunk_df["gc"] == strata]
        rows += strata_df.shape[0]
        strata_df.to_csv(TEMP_DIR + "/" + "{:0.2f}".format(strata),
                        index=False,
                        header=False,
                        mode="a")
    print("chunk_df.shape[0] == rows", chunk_df.shape[0] == rows)
    total_rows += rows

print("------------------ Done splitting ------------------\n")
t1 = time.time()
total = t1 - t0
print("Time:",total,"seconds")
# approx 3 minutes for 1,176,589,961 rows (lib250)


###############################################################################
#### SOME ADJUSTMENTS FOR EMPTY STRATA
# If this .bed has empty strata, we need to accordingly scale
# the number of rows to sample from the remaining non-empty strata.
# 
# note: could still be off by up to 50 (n_strata//2)
# simple example: round(sum([0.5] * 11)) = 6, sum([round(0.5)] * 11) = 11
# however 50 <<< total_rows so this is probably not an issue
# see call to round() in get_sample_n()

t2 = time.time()

print("------------------ Computing sample sizes ------------------")
lines = []
for strata in STRATA:
    strata_file = TEMP_DIR + "/" + "{:0.2f}".format(strata)
    num_lines = sum(1 for line in open(strata_file))
    lines.append(num_lines)
healthy_df["num_lines"] = lines
total_rows = healthy_df["num_lines"].sum() # target total rows
scaling_factor = healthy_df.loc[healthy_df["num_lines"] != 0]["med_frag_fract"].sum()

def get_sample_n(frag_fract, num_lines):
    if frag_fract == 0 or num_lines == 0:
        return 0
    else:
        return round(frag_fract/scaling_factor*total_rows)

healthy_df["sample_n"] = list(map(get_sample_n, healthy_df["med_frag_fract"],healthy_df["num_lines"]))
print("------------------ Done computing sample sizes ------------------")

t3 = time.time()
total = t3 - t2
print("Time:",total,"seconds")
###############################################################################
#### SAMPLE FROM STRATA FILES
# To ensure memory limits when sampling, we limit two values:
# - number of lines to sample at once
# - size of strat_file chunk to read from during sampling

# Because this is evenly drawing from the file w/ replacement,
# we can chunk the sampling, allowing a nice cap to the first value.

# Reading the dataframe in chunks provides a cap to the second value. However
# we want to minimize reading over the file if possible since all this IO gets
# pretty slow

# Two possible strategies:

# STRATEGY 1:
# for each strata:
# 1. get num of lines in file and num of samples
# 2. for each chunk in num of samples
#        3. LINES = sample from all lines in file - range(0,nlines)
#        4. sort LINES
#        5. for each chunk in num of lines in file
#               6. get corresponding lines from LINES
#               7. select corresponding rows from file
#               8. save lines (<tmpdir>/<strata>.sampled)


# STRATEGY 2:
# 1. get num of lines in file and num of samples 
# 2. for each chunk in strata_file:
#        3. compute number of lines to sample from chunk
#        4. for each chunk in num of lines to sample
#               6. sample lines from strata_file chunk
#               7. save lines (<tmpdir>/<strata>.sampled)

# second strategy is easier w/ indexing but acts like a second layer of stratification
# some testing on final statistics will need to be done to see if there was any effect
# TODO: testing of strat 1 vs strat 2

def sample_strata_file(strata_file, target_n, strata):
    num_lines = sum(1 for line in open(strata_file))
    line = 0
    # iterate over strata file in chunks
    while line <= num_lines:
        # compute number of lines to sample from this chunk
        if num_lines - line > CHUNKSIZE:
            sample_range = CHUNKSIZE
        else:
            sample_range = num_lines - line
        sample_n = round(target_n/num_lines*sample_range)
        samples = np.random.choice(sample_range,size=sample_n,replace=True)

        strata_df = pd.read_csv(strata_file,
                        header=None,
                        index_col=False,
                        skiprows=line,
                        nrows=CHUNKSIZE)

        # chunked sampling from strata_file chunk 
        sample_start = 0
        while sample_start < len(samples):
            if sample_start + CHUNKSIZE >= len(samples):
                sample_end = len(samples)
            else:
                sample_end =  sample_start + CHUNKSIZE
            sampled_df = strata_df.iloc[samples[sample_start:sample_end]]
            sampled_df.to_csv(TEMP_DIR + "/{:0.2f}.sampled".format(strata),
                            sep="\t",
                            index=False,
                            header=False,
                            columns=[0,1,2,4,3], #select(chr, start, end, len, gc_strata)
                            mode="a")
            sample_start += CHUNKSIZE
        line += CHUNKSIZE

t4 = time.time()

print("------------------ Sampling strata ------------------")
for strata in STRATA:
    strata_file = TEMP_DIR + "/" + "{:0.2f}".format(strata)
    if not os.path.exists(strata_file):
        print("ERROR?:",strata_file,"not found - skipping!\n")
        continue
    target_n = healthy_df.loc[strata,"sample_n"]
    print("{:0.2f}:".format(strata),target_n)
    # handle empty files
    if target_n == 0:
        open(TEMP_DIR + "/{:0.2f}.sampled".format(strata),"a").close()
    else:
        sample_strata_file(strata_file, target_n, strata)
print("------------------ Done sampling ------------------")

t5 = time.time()
total = t5 - t4
print("Time:",total,"seconds")

###############################################################################
#### MERGE strata files
# https://stackoverflow.com/questions/13613336/how-do-i-concatenate-text-files-in-python

t6 = time.time()

print("------------------ Merging sampled strata + Cleanup ------------------")
sampled_strata = glob.glob(TEMP_DIR + "/*.sampled")

with open(sampled_bed,'wb') as wfd:
    for f in sampled_strata:
        print(f)
        with open(f,'rb') as fd:
            shutil.copyfileobj(fd, wfd)


num_lines = sum(1 for line in open(sampled_bed))
print(num_lines)

#### CLEANUP
try:
    shutil.rmtree(TEMP_DIR)
except OSError as e:
    print("Error: %s - %s." % (e.filename, e.strerror))

print("------------------ Done w/ merging and cleanup ------------------")

t7 = time.time()
total = t7 - t6
print("Time:",total,"seconds")
print("\n\n")

print("Overall time:",t7-t0)


###############################################################################
###############################################################################
# WIP for multithreading and performance notes
###############################################################################
###############################################################################

"""
7,100,413 rows before crash
0.38: 2175208 max strata rows
2,175,208 (~3x) 


lib484 - 9.6G
lib484 - 14G
(1.46x)

lib250 - 6.3G
lib250 - 8.4G
(1.3x)

lib514 - 12G
projected 18G
(18/8 = 2.25x)

not a lot of wiggle room

"""
"""
import multiprocessing as mp

def process_chunk(chunk):
    chunk_df = chunk.round({"gc":2})
    chunk_df["gc_weight"] = chunk_df["gc"].map(healthy_df["med_frag_fract"])
    rows = 0
    # for strata in np.arange(0.0, 1.01, 0.01):
    for strata in np.arange(0.0, 0.02, 0.01):
        strata = round(strata,2) # (!!! ensure floats are properly rounded. will drop rows if removed)
        strata_df = chunk_df.loc[chunk_df["gc"] == strata]
        rows += strata_df.shape[0]
        # strata_chunk_id = "{:0.2f}".format(strata) + "-" + "chunk_" + str(chunk_i)
        # print(chunk_strata_id)
        strata_df.to_csv(TEMP_DIR + "/" + "{:0.2f}".format(strata) + "-" + str(mp.current_process().name),
                        index=False,
                        header=False,)
                        # mode="a")
    # print("chunk_df.shape[0] == rows", chunk_df.shape[0] == rows)
    return rows

print('hi')
reader = pd.read_csv(bed_file,
                        sep="\t", 
                        header=None, 
                        names=COLS, 
                        chunksize=CHUNKSIZE,
                        nrows = CHUNKSIZE*5)

pool = mp.Pool()
funclist = []
for df in reader:
    # process each data frame
    f = pool.apply_async(process_chunk,[df])
    funclist.append(f)

result = 0
for f in funclist:
    result += f.get()

print(result)
print("done")
"""


# import dask.dataframe as dd
# from dask.distributed import Variable

# dask_df = dd.read_csv(bed_file,
#                         sep="\t", 
#                         # usecols=[3], 
#                         header=None, 
#                         names=COLS, 
#                         blocksize=CHUNKSIZE)

# global_var = Variable(name="thread")
# global_var.set(0)

# def process_chunk(chunk_df):
#     chunk_df = chunk.round({"gc":2})
#     chunk_df["gc_weight"] = chunk_df["gc"].map(healthy_df["med_frag_fract"])
#     rows = 0
#     for strata in np.arange(0.0, 1.01, 0.01):
#         strata = round(strata,2) # some float precision nonsense
#         strata_df = chunk_df.loc[chunk_df["gc"] == strata]
#         rows += strata_df.shape[0]
#         thread_id = global_var.get()
#         chunk_strata_id = "{:0.2f}".format(strata) + "-" + "chunk_" + str(thread_id)
#         global_var.set(thread_id+1)
#         print(thread_id)

# result = process_chunk(dask_df)

# print("hi")
# t2 = time.time()
# result = result.compute()
# t3 = time.time()

# print("done")
# print(t3-t2)

# df = pd.DataFrame(columns=["chr","start","end","gc","seq_len","gc_strata"])

# print("hi")
# CHUNKSIZE = 10**6 #(1mb chunks)
# i = 1
# for chunk in pd.read_csv(bed_file,
#                         sep="\t", 
#                         # usecols=[0,1,2,3], 
#                         header=None, 
#                         names=["chr","start","end","gc","seq_len"], 
#                         chunksize=CHUNKSIZE):
#     # if i % 100 == 0:
#     print("chunk number: ", i,flush=True)
#     i += 1
#     print(df.shape)
#     chunk_df = chunk.round({"gc":2})
#     chunk_df["gc_strata"] = chunk_df["gc"].map(healthy_df["med_frag_fract"])
#     print(chunk_df)
#     df = pd.concat([df, chunk_df])

# print(df)

# sampled_df = df.sample(frac=1,replace=True,weights="gc_strata",random_state=1)
# sampled_df.to_csv(sampled_bed,
#     index=False,
#     header=False,
#     sep="\t",
#     columns=["chr","start","end","seq_len","gc"])
