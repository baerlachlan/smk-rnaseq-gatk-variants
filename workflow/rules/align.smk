rule align:
    input:
        unpack(align_inputs),
        idx="resources/genome",
    output:
        aln=temp("results/align/bam/{SAMPLE}_{UNIT}.bam"),
        log="results/align/log/{SAMPLE}_{UNIT}.log",
        log_final="results/align/log/{SAMPLE}_{UNIT}.log.final.out",
    params:
        extra=f"--sjdbOverhang {int(config["read_length"]) - 1} --outSAMtype BAM SortedByCoordinate --twopassMode Basic",
    wrapper:
        "v3.7.0/bio/star/align"


rule align_index:
    input:
        "results/align/bam/{SAMPLE}_{UNIT}.bam",
    output:
        temp("results/align/bam/{SAMPLE}_{UNIT}.bam.bai"),
    wrapper:
        "v3.7.0/bio/samtools/index"


rule align_md5:
    input:
        align_md5_inputs(),
    output:
        "results/align/bam/md5.txt",
    shell:
        """
        md5sum {input} > {output}
        """
