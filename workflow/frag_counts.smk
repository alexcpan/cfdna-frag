rule gc_distro:
    input:
        frag = config["data_dir"] + "/frag/{library_id}_frag.bed",
    params:
        config["r_lib_loads"],
    output:
        config["data_dir"] + "/frag/{library_id}_gc_distro.csv"
    script:
        "scripts/gc_distro.R"

rule make_healthy_gc_summary:
    output:
        healthy_med = config["data_dir"] + "/frag/healthy_med.rds"
    script:
        "scripts/make_healthy_gc_summary.R"

rule sample_frags_by_gc:
    input:
        healthy_med = config["data_dir"] + "/frag/healthy_med.rds",
        frag_bed = config["data_dir"] + "/frag/{library_id}_frag.bed"
    output:
        config["data_dir"] + "/frag/{library_id}_norm_frag.bed"
    script:
        "scripts/sample_frags_by_gc.R"

rule frag_window_sum:
    input:
        frag = config["data_dir"] + "/frag/{library_id}_norm_frag.bed",
    output:
        short = config["data_dir"] + "/frag/{library_id}_norm_short.bed",
        long = config["data_dir"] + "/frag/{library_id}_norm_long.bed",
    shell:
        """
        workflow/scripts/frag_window_sum.sh {input.frag} \
                                            {output.short} \
                                            {output.long}
        """

rule frag_window_int:
    input:
        short = config["data_dir"] + "/frag/{library_id}_norm_short.bed",
        long = config["data_dir"] + "/frag/{library_id}_norm_long.bed",
        matbed = config["data_dir"] + "/ref/mathios_chrom_bins.bed",
    output:
        cnt_long_tmp = config["data_dir"] + "/frag/{library_id}_cnt_long.tmp",
        cnt_short_tmp = config["data_dir"] + "/frag/{library_id}_cnt_short.tmp",
        cnt_long = config["data_dir"] + "/frag/{library_id}_cnt_long.bed",
        cnt_short = config["data_dir"] + "/frag/{library_id}_cnt_short.bed",
    shell:
        """
        bedtools intersect -c -a {input.matbed} -b {input.long} > {output.cnt_long_tmp}
        awk '{{print FILENAME (NF?"\t":"") $0}}' {output.cnt_long_tmp} |
        sed 's/^.*lib/lib/g' |
        sed 's/_cnt_/\t/g' |
        sed 's/.tmp//g' |
        awk 'BEGIN {{OFS="\t"}}; {{print $1,$2,$3,$4,$5,$10}}' > {output.cnt_long}
        bedtools intersect -c -a {input.matbed} -b {input.short} > {output.cnt_short_tmp}
        awk '{{print FILENAME (NF?"\t":"") $0}}' {output.cnt_short_tmp} |
        sed 's/^.*lib/lib/g' |
        sed 's/_cnt_/\t/g' |
        sed 's/.tmp//g' |
        awk 'BEGIN {{OFS="\t"}}; {{print $1,$2,$3,$4,$5,$10}}' > {output.cnt_short}
        """

rule count_merge:
    input:
        expand(config["data_dir"] + "/frag/{library_id}_cnt_{length}.bed", library_id=ALLLIB, length=["short", "long"])
    output:
        config["data_dir"] + "/frag/frag_counts.tsv"
    shell:
        """
        cat {input} > {output}
        """

rule count_scale:
    input:
    output:
    script:
        "scripts/count_scale.R"
