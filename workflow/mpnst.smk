# cfdna_wgs_frag_bam_inputs = config["dir"]["data"]["cfdna_wgs"] + "/bam/raw"
# cfdna_wgs_frag_filt_bams  = config["dir"]["data"]["cfdna_wgs"] + "/bam/frag"
# cfdna_wgs_frag_bam_inputs = config["dir"]["data"]["cfdna_wgs"] + "/bam"
# cfdna_wgs_frag_filt_bams  = config["dir"]["data"]["cfdna_wgs"] + "/bam"
cfdna_wgs_frag_bam_inputs = "/data/panalc/MPNST/cfdna-wgs-bams"
cfdna_wgs_frag_filt_bams  = "/data/panalc/MPNST/cfdna-wgs-bams"
cfdna_wgs_frag_filt_sort_bams  = config["dir"]["data"]["cfdna_wgs"] + "/bam"
cfdna_wgs_frag_beds =       config["dir"]["data"]["cfdna_wgs"] + "/frag"
cfdna_wgs_distros =         config["dir"]["data"]["cfdna_wgs"] + "/distro"
cfdna_wgs_logs =            config["dir"]["data"]["cfdna_wgs"] + "/logs"
cfdna_wgs_frag_len =        config["dir"]["data"]["cfdna_wgs"] + "/len"
cfdna_wgs_frag_cnt =        config["dir"]["data"]["cfdna_wgs"] + "/bed-frag-cnt"


print(cfdna_wgs_frag_beds)

import pandas as pd
df = pd.read_csv("/data/panalc/cfdna-frag/config/new_library.csv")
dx = ["healthy","plexiform","an","annubp","mpnst"]
df = df.loc[
    (df["seq"] == "wgs") &
    (df["isolation"] == "cfdna") &
    (df["current_dx"].isin(dx)) &
    (df["processed"] == 1)
]
print(df["current_dx"].value_counts())
LIBRARIES = df["library"].to_list()
LIBRARIES_HEALTHY = df.loc[df["current_dx"] == "healthy"]["library"].to_list()

rule all:
    input:
        expand(cfdna_wgs_frag_filt_bams + "/{library}_filt.bam", library = LIBRARIES),
        expand(cfdna_wgs_frag_beds      + "/{library_id}_frag.bed", library_id = LIBRARIES),
        expand(cfdna_wgs_distros        + "/{library_id}_gc_distro.csv", library_id = LIBRARIES),
        # cfdna_wgs_distros + "/healthy_med.rds",
        cfdna_wgs_distros + "/healthy_med.csv",
        expand(cfdna_wgs_frag_beds      + "/{library_id}_frag_sampled.bed", library_id = LIBRARIES),
        expand(cfdna_wgs_frag_len      + "/{library_id}_norm_short.bed", library_id = LIBRARIES),
        expand(cfdna_wgs_frag_cnt + "/{library_id}_cnt_long.bed", library_id = LIBRARIES),

include: config["dir"]["repo"]["cfdna_wgs"] + "/workflow/frag_bed.smk"
