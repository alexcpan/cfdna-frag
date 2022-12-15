#! /bin/bash
set -o errexit   # abort on nonzero exitstatus
set -o nounset   # abort on unbound variable
set -o pipefail  # don't hide errors within pipes

input_bam="${1}"
fasta="${2}"
frag_dir="${3}"
lib_id="${4}"

tmp_dir="${frag_dir}/${lib_id}_tmp"
tmp_bed="${tmp_dir}/${lib_id}_tmp.bed"
out_bed="${frag_dir}/${lib_id}_frag.bed"

mkdir -p $tmp_dir

# convert to bed
bedtools bamtobed -bedpe -i $input_bam |
awk '$1==$4 {print $0}' | awk '$2 < $6 {print $0}' |
awk -v OFS='\t' '{print $1,$2,$6}' > \
$tmp_bed

wc -l $tmp_bed

# split bed into 10,000 line chunks
split -l 10000 -a 4 $tmp_bed \
                    "${tmp_bed}.pe."

pe_files="${tmp_bed}.pe.*"
ls -l $pe_files | wc -l

# use pyBedTools to calculate nucleotide content over each chunk
python /data/panalc/cfdna-frag/workflow/scripts/readtofrag.py \
--bed_dir=$tmp_dir \
--fasta=$fasta

nc_files="${tmp_dir}/*.nc"
ls -l $nc_files | wc -l

# merge chunks back tgoether
awk 'FNR>1' $nc_files |
awk -v OFS='\t' '{print $1,$2,$3,$5,$12}' > \
$out_bed

wc -l $out_bed

rm -rf $tmp_dir
