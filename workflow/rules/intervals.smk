rule gtf_to_bed:
    input:
        gtf="resources/annotation.gtf",
    output:
        bed=temp("resources/exons.bed"),
    conda:
        "../envs/r-base.yml"
    script:
        "../scripts/gtf_to_bed.R"


rule bed_to_intervals:
    input:
        bed="resources/exons.bed",
        dict="resources/genome.dict",
    output:
        temp("resources/exons.interval_list"),
    wrapper:
        "v3.7.0/bio/picard/bedtointervallist"
