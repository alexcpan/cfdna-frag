# For each library, makes a csv with columns of library_id, gc_strata, and fract_frags
rule gc_distro:
    container:
        config["container"]["cfdna_wgs"],
    input:
        cfdna_wgs_frag_beds + "/{library_id}_frag.bed",
    log:
        cfdna_wgs_logs + "/{library_id}_gc_distro.log",
    output:
        cfdna_wgs_distros + "/{library_id}_gc_distro.csv"
    params:
        script = config["dir"]["scripts"]["cfdna_wgs"] + "/gc_distro.R",
    shell:
        """
        Rscript {params.script} \
        {input} \
        {output} \
        > {log} 2>&1
        """

# Make tibble of gc_strata and median fraction of fragments from healthy samples
rule make_healthy_gc_summary:
    container:
        config["container"]["cfdna_wgs"],
    log:
        cfdna_wgs_logs + "/make_healthy_gc_summary.log",
    output:
        cfdna_wgs_distros + "/healthy_med.rds"
    params:
        distro_dir = cfdna_wgs_distros,
        healthy_libs_str = LIBRARIES_HEALTHY,
        script = config["dir"]["scripts"]["cfdna_wgs"] + "/make_healthy_gc_summary.R",
    shell:
        """
        Rscript {params.script} \
        {params.distro_dir} \
        "{params.healthy_libs_str}" \
        {output} > {log} 2>&1
        """

rule sample_frags_by_gc:
    container:
        config["container"]["cfdna_wgs"],
    input:
        healthy_med = cfdna_wgs_distros + "/healthy_med.rds",
        frag_bed = cfdna_wgs_frag_beds + "/{library_id}_frag.bed",
    log:
        cfdna_wgs_logs + "/{library_id}_sample_frags_by_gc.log",
    output:
        cfdna_wgs_frag_beds + "/{library_id}_frag_sampled.bed",
    params:
        script = config["dir"]["scripts"]["cfdna_wgs"] + "/sample_frags_by_gc.R",
    shell:
        """
        Rscript {params.script} \
        {input.healthy_med} \
        {input.frag_bed} \
        {output} > {log} 2>&1
        """

rule frag_window_sum:
    container:
        config["container"]["cfdna_wgs"],
    input:
        cfdna_wgs_frag_beds + "/{library_id}_norm_frag.bed",
    log:
        cfdna_wgs_logs + "/{library_id}_frag_window_sum.log",
    output:
        short = cfdna_wgs_frag_len + "/{library_id}_norm_short.bed",
        long = cfdna_wgs_frag_len + "/{library_id}_norm_long.bed",
    params:
        script = config["dir"]["scripts"]["cfdna_wgs"] + "/frag_window_sum.sh",
    shell:
        """
        {params.script} \
        {input.frag} \
        {output.short} \
        {output.long} &> {log}
        """
