rule filter_alignments:
    input:
        bam = cfdna_wgs_frag_bam_inputs + "/{library}.bam",
        keep_bed = config["files"]["cfdna_wgs_frag_keep_bed"],
    params:
        script = config["dir"]["scripts"]["cfdna_wgs"] + "/filter_alignments.sh",
        temp_dir = config["dir"]["data"]["cfdna_wgs"] + "/tmp",
        threads = config["threads"]["cfdna_wgs"],
    resources:
        mem_mb=5000
    output:
        cfdna_wgs_frag_filt_bams + "/{library}_filt.bam",
    container:
        config["container"]["cfdna_wgs"],
    shell:
        """
        {params.script} \
        {input.bam} \
        {input.keep_bed} \
        {params.temp_dir} \
        {params.threads} \
        {output}
        """

# new
# technically we start from here since the filtering step was performed in the wgs-cfdna pipline
rule sort_bam:
    input:
        bam = cfdna_wgs_frag_filt_bams + "/{library}_filt.bam",
    params:
        script = config["dir"]["scripts"]["cfdna_wgs"] + "/sort_bam.sh",
        temp_dir = config["dir"]["data"]["cfdna_wgs"] + "/tmp",
    resources:
        time   = "1:00:00",
        mem_gb = "10g",
        cpus   = "20",
    output:
        cfdna_wgs_frag_filt_sort_bams + "/{library}_filt.sort.bam",
    container:
        config["container"]["cfdna_wgs"],
    shell:
        """
        {params.script} \
        {input.bam} \
        {params.temp_dir} \
        {resources.cpus} \
        {output}
        """

## modified
# 1. bam is converted to bed
# 2. bed gets split into 10,000 line chunks
# 3. pyBedTools calculates nucleotide content over each chunk in parallel
# (automatically uses optimal number of threads)
# 4. merges nuc content files into a single bed
rule read_to_frag_bed:
    input:
        # cfdna_wgs_frag_filt_bams + "/{library_id}_filt.bam",
        cfdna_wgs_frag_filt_sort_bams + "/{library_id}_filt.sort.bam",
    params:
        # fasta = config["dir"]["data"]["cfdna_wgs"] + "/inputs/chr19.fa",
        # script = config["dir"]["scripts"]["cfdna_wgs"] + "/read_to_frag_bed.sh",
        fasta    = "/data/panalc/cfdna-frag/MPNST_frag/inputs/hg38.fa",
        script   = config["dir"]["scripts"]["cfdna_wgs"] + "/readtofrag.sh",
        frag_dir = cfdna_wgs_frag_beds,
        lib_id   = "{library_id}",
    output:
        cfdna_wgs_frag_beds + "/{library_id}_frag.bed",
    resources:
        time   = "1:00:00",
        mem_gb = "10g",
        cpus   = "20",
    # container:
    #     config["container"]["cfdna_wgs"]
    container: None
    shell:
        # {params.script} \
	    # {input} \
        # {params.fasta} \
        # {output}
        """
        {params.script} \
	    {input} \
        {params.fasta} \
        {params.frag_dir} \
        {params.lib_id}
        """

## modified
# -> instead of reading over a whole file, compute sums over chunks and then average at the end
# For each library, makes a csv with columns of library_id, gc_strata, and fract_frags
rule gc_distro:
    # container:
    #     config["container"]["cfdna_wgs"],
    container: None
    input:
        cfdna_wgs_frag_beds + "/{library_id}_frag.bed",
    log:
        cfdna_wgs_logs + "/{library_id}_gc_distro.log",
    output:
        cfdna_wgs_distros + "/{library_id}_gc_distro.csv"
    params:
        # script = config["dir"]["scripts"]["cfdna_wgs"] + "/gc_distro.R",
        script = config["dir"]["scripts"]["cfdna_wgs"] + "/gc_distro.py",
    shell:
        # """
        # Rscript {params.script} \
        # {input} \
        # {output} \
        # > {log} 2>&1
        # """
        """
        python3 {params.script} \
        --bed_file {input} \
        --distro_file {output} \
        > {log} 2>&1
        """

## modified
# -> save to a .csv instead of a .rds so that we can read it from python later
# Make tibble of gc_strata and median fraction of fragments from healthy samples

### XXX: this gives a distribution that doesn't add to 1 since it's based on medians.
# we may want to normalize it to 1 instead? not really a problem though

rule make_healthy_gc_summary:
    container:
        config["container"]["cfdna_wgs"],
    input:
        expand(cfdna_wgs_distros + "/{library_id}_gc_distro.csv", library_id = LIBRARIES),
    log:
        cfdna_wgs_logs + "/make_healthy_gc_summary.log",
    output:
        # cfdna_wgs_distros + "/healthy_med.rds"
        cfdna_wgs_distros + "/healthy_med.csv"
    params:
        distro_dir = cfdna_wgs_distros,
        healthy_libs_str = LIBRARIES_HEALTHY,
        # script = config["dir"]["scripts"]["cfdna_wgs"] + "/make_healthy_gc_summary.R",
        script = config["dir"]["scripts"]["cfdna_wgs"] + "/make_healthy_gc_summary_csv.R",
    shell:
        """
        Rscript {params.script} \
        {params.distro_dir} \
        "{params.healthy_libs_str}" \
        {output} \
        > {log} 2>&1
        """

#  modified
# previously samples over the entire bed file which was too big to read into memory on biowulf
# 1. read a bam in chunks, storing each fragment in the chunk to a respective strata file
# done in parallel, so each thread creates its own set of strata files
# theoretically creates (n_strata x n_chunks) files where n_strata = 101 (0.00 to 1.00)
# 2. merges the chunks for each strata together to get 101 strata files
# 3. compute # of fragments to sample (w/ replacement) from each chunk
# 4. sample w/ replacement from each strata file -> get 101 sampled files
# 5. merge sampled strata files into a single file
# (6. sort merged file by genomic coordinates) TODO
# step 6 requires external merging which is a pain to implement
# for whatever reason, using unix sort function still OOMs even with temp dir option

rule sample_frags_by_gc:
    container:
        config["container"]["cfdna_wgs"],
    input:
        # healthy_med = cfdna_wgs_distros + "/healthy_med.rds",
        healthy_med = cfdna_wgs_distros + "/healthy_med.csv",
        frag_bed = cfdna_wgs_frag_beds + "/{library_id}_frag.bed",
    log:
        cfdna_wgs_logs + "/{library_id}_sample_frags_by_gc.log",
    output:
        cfdna_wgs_frag_beds + "/{library_id}_frag_sampled.bed",
    params:
        # script = config["dir"]["scripts"]["cfdna_wgs"] + "/sample_frags_by_gc.R",
        script = config["dir"]["scripts"]["cfdna_wgs"] + "/sample_frags_by_gc.py",
    shell:
        # """
        # Rscript {params.script} \
        # {input.healthy_med} \
        # {input.frag_bed} \
        # {output} > {log} 2>&1
        # """
        """
        python3 {params.script} \
        --healthy_csv {input.healthy_med} \
        --bed_file {input.frag_bed} \
        --sampled_bed {output} > {log} 2>&1
        """


rule frag_window_sum:
    container:
        config["container"]["cfdna_wgs"],
    input:
        cfdna_wgs_frag_beds + "/{library_id}_frag_sampled.bed",
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
        {input} \
        {output.short} \
        {output.long} &> {log}
        """

# modified
rule frag_window_int:
    input:
        short = cfdna_wgs_frag_len + "/{library_id}_norm_short.bed",
        long = cfdna_wgs_frag_len + "/{library_id}_norm_long.bed",
        matbed = config["files"]["cfdna_wgs_frag_keep_bed"]
        # matbed = "test/inputs/keep.bed",
    output:
        cnt_long_tmp = temp(cfdna_wgs_frag_cnt + "/{library_id}_cnt_long.tmp"),
        cnt_short_tmp = temp(cfdna_wgs_frag_cnt + "/{library_id}_cnt_short.tmp"),
        cnt_long = cfdna_wgs_frag_cnt + "/{library_id}_cnt_long.bed",
        cnt_short = cfdna_wgs_frag_cnt + "/{library_id}_cnt_short.bed",
    resources:
        time   = "1:00:00",
        mem_gb = "30g",
        cpus   = "2",
    shell:
        """
        module load bedtools
        bedtools intersect -c -a {input.matbed} -b {input.long} > {output.cnt_long_tmp}
        awk '{{print FILENAME (NF?"\t":"") $0}}' {output.cnt_long_tmp} |
        sed 's/^.*lib/lib/g' |
        sed 's/_.*tmp//g' > {output.cnt_long}
        bedtools intersect -c -a {input.matbed} -b {input.short} > {output.cnt_short_tmp}
        awk '{{print FILENAME (NF?"\t":"") $0}}' {output.cnt_short_tmp} |
        sed 's/^.*lib/lib/g' |
        sed 's/_.*tmp//g' > {output.cnt_short}
        """
