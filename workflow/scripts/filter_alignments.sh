# Function
filter_bams(){
    # Filter to mapq 30 and limit to keep.bed genomic regions
    samtools view -@ $1 -b -h -L $2 -o - -q 30 $3 |
    samtools sort -@ $1 -n -o $4 -T $5 -
    }

# Snakemake variables
input_in_bam="$1"
input_keep_bed="$2"
params_temp_dir="$3"
params_threads="$4"
output_filt_bam="$5"

# Run command
filter_bams "$params_threads" "$input_keep_bed" "$input_in_bam" "$output_filt_bam" $params_temp_dir
