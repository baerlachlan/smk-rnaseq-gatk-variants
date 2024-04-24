rule assign_read_groups:
    input:
        "results/align/bam/{SAMPLE}_{UNIT}.bam",
    output:
        temp("results/assign_read_groups/bam/{SAMPLE}_{UNIT}.bam"),
    params:
        extra=lambda w: f"--RGPU {w.SAMPLE}_{w.UNIT} --RGSM {w.SAMPLE} --RGPL ILLUMINA --RGLB null --SORT_ORDER coordinate",
    wrapper:
        "v3.7.0/bio/picard/addorreplacereadgroups"

rule assign_read_groups_index:
    input:
        "results/assign_read_groups/bam/{SAMPLE}_{UNIT}.bam",
    output:
        temp(temp("results/assign_read_groups/bam/{SAMPLE}_{UNIT}.bam.bai")),
    wrapper:
        "v3.7.0/bio/samtools/index"

rule assign_read_groups_md5:
    input:
        assign_read_groups_md5_inputs(),
    output:
        "results/assign_read_groups/bam/md5.txt",
    shell:
        """
        md5sum {input} > {output}
        """