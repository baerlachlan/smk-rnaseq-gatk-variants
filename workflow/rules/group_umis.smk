rule group_umis_se:
    input:
        bam="results/assign_read_groups/bam/{SAMPLE}_{UNIT}.bam",
        bam_idx="results/assign_read_groups/bam/{SAMPLE}_{UNIT}.bam.bai",
    output:
        bam=temp("results/group_umis/bam/{SAMPLE}_{UNIT}_single_end.bam"),
        bam_idx=temp("results/group_umis/bam/{SAMPLE}_{UNIT}_single_end.bam.bai"),
    conda:
        "../envs/umi-tools.yml"
    shell:
        """
        umi_tools group \
            -I {input.bam} \
            -S {output.bam} \
            --temp-dir=$TMPDIR \
            --output-bam \
            --method=unique \
            --extract-umi-method=read_id \
            --umi-separator=":"

        samtools index {output.bam}
        """


rule group_umis_pe:
    input:
        bam="results/assign_read_groups/bam/{SAMPLE}_{UNIT}.bam",
        bam_idx="results/assign_read_groups/bam/{SAMPLE}_{UNIT}.bam.bai",
    output:
        bam=temp("results/group_umis/bam/{SAMPLE}_{UNIT}_paired_end.bam"),
        bam_idx=temp("results/group_umis/bam/{SAMPLE}_{UNIT}_paired_end.bam.bai"),
    conda:
        "../envs/umi-tools.yml"
    shell:
        """
        umi_tools group \
            -I {input.bam} \
            -S {output.bam} \
            --temp-dir=$TMPDIR \
            --output-bam \
            --method=unique \
            --extract-umi-method=read_id \
            --umi-separator=":" \
            --paired

        samtools index {output.bam}
        """


rule group_umis_md5:
    input:
        group_umis_md5_inputs(),
    output:
        "results/group_umis/bam/md5.txt",
    conda:
        "../envs/parallel.yml"
    shell:
        """
        echo "{input}" | tr " " "\n" | parallel -j {threads} md5sum  > {output}
        """