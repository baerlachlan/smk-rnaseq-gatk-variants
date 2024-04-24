rule mark_duplicates:
    input:
        unpack(mark_duplicates_inputs),
    output:
        bam=temp("results/mark_duplicates/bam/{SAMPLE}.bam"),
        metrics="results/mark_duplicates/log/{SAMPLE}.metrics.txt",
    wrapper:
        "v3.7.0/bio/picard/markduplicates"

rule mark_duplicates_umi:
    input:
        unpack(mark_duplicates_inputs),
    output:
        bam=temp("results/mark_duplicates/bam/{SAMPLE}_umi.bam"),
        metrics="results/mark_duplicates/log/{SAMPLE}_umi.metrics.txt",
    params:
        extra="--BARCODE_TAG BX"
    wrapper:
        "v3.7.0/bio/picard/markduplicates"

rule mark_duplicates_md5:
    input:
        mark_duplicates_md5_inputs(),
    output:
        "results/mark_duplicates/bam/md5.txt",
    shell:
        """
        md5sum {input} > {output}
        """