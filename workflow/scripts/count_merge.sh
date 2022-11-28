array=($1)
output=$2

if [ -f $output ]; then \rm $output; fi

for file in ${array[@]}; do
    awk '{{print FILENAME (NF?"\t":"") $0}}' $file |
        sed 's/^.*lib/lib/g' |
        sed 's/_.*_/\t/g' |
        sed 's/\.bed//g' >> $output
done
