rule bqsr_firstpass:
    input:
        bam="results/split_n_cigar_reads/bam/{SAMPLE}.bam",
        ref="resources/genome.fa",
        ref_idx="resources/genome.fa.fai",
        dict="resources/genome.dict",
        known="resources/known_variants.vcf.gz",
        known_idx="resources/known_variants.vcf.gz.tbi",
    output:
        recal_table=temp("results/bqsr/recal/{SAMPLE}.grp"),
    wrapper:
        "v4.0.0/bio/gatk/baserecalibrator"


rule bqsr_apply:
    input:
        bam="results/split_n_cigar_reads/bam/{SAMPLE}.bam",
        ref="resources/genome.fa",
        ref_idx="resources/genome.fa.fai",
        dict="resources/genome.dict",
        recal_table="results/bqsr/recal/{SAMPLE}.grp",
    output:
        bam=temp("results/bqsr/bam/{SAMPLE}.bam"),
        bam_idx=temp("results/bqsr/bam/{SAMPLE}.bai"),
    params:
        extra="--add-output-sam-program-record",
    wrapper:
        "v4.0.0/bio/gatk/applybqsr"


rule bqsr_apply_md5:
    input:
        expand("results/bqsr/bam/{SAMPLE}.bam", SAMPLE=samples["sample"]),
    output:
        "results/bqsr/bam/md5.txt",
    conda:
        "../envs/parallel.yml"
    shell:
        """
        echo "{input}" | tr " " "\n" | parallel -j {threads} md5sum  > {output}
        """
