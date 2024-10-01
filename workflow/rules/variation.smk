rule variation_call:
    input:
        bam="results/bqsr/bam/{SAMPLE}.bam",
        bam_idx="results/bqsr/bam/{SAMPLE}.bai",
        ref="resources/genome.fa",
        ref_idx="resources/genome.fa.fai",
        dict="resources/genome.dict",
        intervals="resources/exons.interval_list",
        known="resources/known_variants.vcf.gz",
        known_idx="resources/known_variants.vcf.gz.tbi",
    output:
        gvcf=temp("results/variation/gvcf/{SAMPLE}.g.vcf"),
        gvcf_idx=temp("results/variation/gvcf/{SAMPLE}.g.vcf.idx"),
    params:
        extra="-dont-use-soft-clipped-bases --standard-min-confidence-threshold-for-calling 20",
    wrapper:
        "v4.0.0/bio/gatk/haplotypecaller"


rule variation_combine:
    input:
        gvcfs=expand("results/variation/gvcf/{SAMPLE}.g.vcf", SAMPLE=samples["sample"]),
        ref="resources/genome.fa",
        ref_idx="resources/genome.fa.fai",
        dict="resources/genome.dict",
    output:
        gvcf=temp("results/variation/gvcf/all.g.vcf"),
        gvcf_idx=temp("results/variation/gvcf/all.g.vcf.idx"),
    wrapper:
        "v4.0.0/bio/gatk/combinegvcfs"


rule variation_genotype:
    input:
        gvcf="results/variation/gvcf/all.g.vcf",
        ref="resources/genome.fa",
        ref_idx="resources/genome.fa.fai",
        dict="resources/genome.dict",
    output:
        vcf=temp("results/variation/vcf/all.vcf"),
        vcf_idx=temp("results/variation/vcf/all.vcf.idx"),
    wrapper:
        "v4.0.0/bio/gatk/genotypegvcfs"


rule variation_select:
    input:
        vcf="results/variation/vcf/all.vcf",
        ref="resources/genome.fa",
        ref_idx="resources/genome.fa.fai",
        dict="resources/genome.dict",
    output:
        vcf=temp("results/variation/vcf/snvs.vcf"),
        vcf_idx=temp("results/variation/vcf/snvs.vcf.idx"),
    params:
        extra="--select-type-to-include SNP",
    wrapper:
        "v4.0.0/bio/gatk/selectvariants"


rule variation_filter:
    input:
        vcf="results/variation/vcf/snvs.vcf",
        ref="resources/genome.fa",
        ref_idx="resources/genome.fa.fai",
        dict="resources/genome.dict",
    output:
        vcf="results/variation/vcf/snvs.filtered.vcf",
        vcf_idx="results/variation/vcf/snvs.filtered.vcf.idx",
    params:
        filters={
            "FS": "FS > 60.0",
            "QD": "QD < 2.0",
            "MQ": "MQ < 40.0",
            "SOR": "SOR > 4.0",
            "MQRankSum": "MQRankSum < -12.5",
            "ReadPosRankSum": "ReadPosRankSum < -8.0",
        },
    wrapper:
        "v4.0.0/bio/gatk/variantfiltration"


rule variation_cleanup:
    input:
        vcf="results/variation/vcf/snvs.filtered.vcf",
        ref="resources/genome.fa",
        ref_idx="resources/genome.fa.fai",
        dict="resources/genome.dict",
    output:
        vcf="results/variation/vcf/snvs.final.vcf",
        vcf_idx="results/variation/vcf/snvs.final.vcf.idx",
    params:
        extra="--exclude-filtered",
    wrapper:
        "v4.0.0/bio/gatk/selectvariants"


rule variation_gvcf_md5:
    input:
        singles=expand(
            "results/variation/gvcf/{SAMPLE}.g.vcf", SAMPLE=samples["sample"]
        ),
        combined="results/variation/gvcf/all.g.vcf",
    output:
        "results/variation/gvcf/md5.txt",
    conda:
        "../envs/parallel.yml"
    shell:
        """
        echo "{input}" | tr " " "\n" | parallel -j {threads} md5sum  > {output}
        """


rule variation_vcf_md5:
    input:
        [
            "results/variation/vcf/all.vcf",
            "results/variation/vcf/snvs.vcf",
            "results/variation/vcf/snvs.filtered.vcf",
            "results/variation/vcf/snvs.final.vcf",
        ],
    output:
        "results/variation/vcf/md5.txt",
    conda:
        "../envs/parallel.yml"
    shell:
        """
        echo "{input}" | tr " " "\n" | parallel -j {threads} md5sum  > {output}
        """
