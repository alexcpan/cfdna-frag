cfdna_wgs_frag_bam_inputs = config["dir"]["data"]["cfdna_wgs"] + "/bam/raw"
cfdna_wgs_frag_filt_bams  = config["dir"]["data"]["cfdna_wgs"] + "/bam/frag"
cfdna_wgs_frag_beds =       config["dir"]["data"]["cfdna_wgs"] + "/frag"
cfdna_wgs_distros =         config["dir"]["data"]["cfdna_wgs"] + "/distro"
cfdna_wgs_logs =            config["dir"]["data"]["cfdna_wgs"] + "/logs"
cfdna_wgs_frag_len =        config["dir"]["data"]["cfdna_wgs"] + "/len"
cfdna_wgs_frag_cnt =        config["dir"]["data"]["cfdna_wgs"] + "/bed-frag-cnt"

LIBRARIES = ["lib001", "lib002"]

LIBRARIES_HEALTHY = ["lib001", "lib002"]

rule all:
    input:
        expand(cfdna_wgs_frag_filt_bams + "/{library}_filt.bam", library = LIBRARIES),
        expand(cfdna_wgs_frag_beds      + "/{library_id}_frag.bed", library_id = LIBRARIES),
        expand(cfdna_wgs_distros        + "/{library_id}_gc_distro.csv", library_id = LIBRARIES),
        cfdna_wgs_distros + "/healthy_med.rds",
        expand(cfdna_wgs_frag_beds      + "/{library_id}_frag_sampled.bed", library_id = LIBRARIES),
        expand(cfdna_wgs_frag_len      + "/{library_id}_norm_short.bed", library_id = LIBRARIES),
        expand(cfdna_wgs_frag_cnt + "/{library_id}_cnt_long.bed", library_id = LIBRARIES),

include: config["dir"]["repo"]["cfdna_wgs"] + "/workflow/frag_bed.smk"
