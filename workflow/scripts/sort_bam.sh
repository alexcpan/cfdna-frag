#!/usr/bin/env bash

#### SEE filter_alignments.sh ####
# Function
# sort_bam(){
#     # # Filter to mapq 30 and limit to keep.bed genomic regions
#     # samtools view -@ $1 -b -h -L $2 -o - -q 30 $3 |
#     # samtools sort -@ $1 -n -o - -T $5 - |
#     # samtools fixmate -@ $1 - - |
#     samtools sort -@ $1 -n -o $4 -T $5 $3
#     }

# # Snakemake variables
# input_in_bam="$1"
# input_keep_bed="$2"
# params_temp_dir="$3"
# params_threads="$4"
# output_filt_bam="$5"

# # Run command
# sort_bam "$params_threads" "$input_keep_bed" "$input_in_bam" "$output_filt_bam" $params_temp_dir

####################################################################################

input_bam="${1}"
tmp_dir="${2}"
threads="${3}"
out_bam="${4}"

sambamba sort -t $threads \
    --tmpdir=$tmp_dir \
    --show-progress \
    --sort-by-name \
    --out=/dev/stdout \
    $input_bam |
samtools fixmate -@ $threads - - |
sambamba sort -t $threads \
    --tmpdir=$tmp_dir \
    --show-progress \
    --sort-by-name \
    --out=$out_bam \
    /dev/stdin
