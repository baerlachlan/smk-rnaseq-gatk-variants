rule split_n_cigar_reads:
    input:
        unpack(split_n_cigar_reads_inputs),
        ref="resources/genome.fa",
        ref_idx="resources/genome.fa.fai",
    output:
        temp("results/split_n_cigar_reads/bam/{SAMPLE}.bam"),
        bam_idx=temp("results/split_n_cigar_reads/bam/{SAMPLE}.bai"),
    wrapper:
        "v4.0.0/bio/gatk/splitncigarreads"


rule split_n_cigar_reads_md5:
    input:
        expand("results/split_n_cigar_reads/bam/{SAMPLE}.bam", SAMPLE=samples["sample"]),
    output:
        "results/split_n_cigar_reads/bam/md5.txt",
    conda:
        "../envs/parallel.yml"
    shell:
        """
        echo "{input}" | tr " " "\n" | parallel -j {threads} md5sum  > {output}
        """
