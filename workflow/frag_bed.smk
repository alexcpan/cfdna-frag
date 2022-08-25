rule filter_alignments:
    input:
        bam = cfdna_wgs_frag_bam_inputs + "/{library_id}.bam",
        keep_bed = config["files"]["cfdna_wgs_frag_keep_bed"],
    params:
        temp_dir = config["dir"]["data"]["cfdna_wgs"] + "/tmp",
        threads = config["threads"]["cfdna_wgs"],
    resources:
        mem_mb=5000
    output:
        cfdna_wgs_frag_filt_bams + "/{library_id}_filt.bam",
    container:
        config["container"]["cfdna_wgs"],
    shell:
        """
        workflow/scripts/filter_alignments.sh {input.bam} \
                                              {input.keep_bed} \
                                              {params.temp_dir} \
                                              {params.threads} \
                                              {output}
        """

rule read_to_frag_bed:
    input:
        cfdna_wgs_frag_filt_bams + "/{library_id}_filt.bam",
    params:
        fasta = config["files"]["cfdna_wgs_genome_fasta"],
    output:
        cfdna_wgs_frag_beds + "/{library_id}_frag.bed",
    resources:
        mem_mb=5000
    container:
        config["container"]["cfdna_wgs"]
    shell:
        """
        workflow/scripts/read_to_frag_bed.sh {input} \
                                             {params.fasta} \
                                             {output}
        """
