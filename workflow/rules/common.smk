import pandas as pd
from snakemake.utils import min_version, validate

min_version("8.4.2")


configfile: "config/config.yaml"


samples = pd.read_csv(config["samples"], sep="\t", dtype={"sample": str}).set_index(
    "sample", drop=False
)
validate(samples, "../schemas/samples.schema.yml")

units = pd.read_csv(
    config["units"], sep="\t", dtype={"sample": str, "unit": str}
).set_index(["sample", "unit"], drop=False)
validate(units, "../schemas/units.schema.yml")

####
## Helper functions
####


def is_paired_end(sample):
    sample_units = units.loc[sample]
    fq2_null = sample_units["fq2"].isnull()
    paired = ~fq2_null
    all_paired = paired.all()
    all_single = (~paired).all()
    assert (
        all_single or all_paired
    ), f"all units for sample {sample} must be single or paired end"
    return all_paired


if all(is_paired_end(i) for i in samples["sample"]):
    pair_tags = ["R1", "R2"]
else:
    pair_tags = ["R0"]

####
## Wildcard constraints
####


wildcard_constraints:
    SAMPLE="|".join(samples["sample"]),
    UNIT="|".join(units["unit"]),
    PAIRTAG="|".join(pair_tags),


####
## Input functions
####


def fastqc_raw_inputs(wildcards):
    unit = units.loc[wildcards.SAMPLE, wildcards.UNIT]
    if is_paired_end(wildcards.SAMPLE):
        if wildcards.PAIRTAG == pair_tags[0]:
            return f"{unit.fq1}"
        elif wildcards.PAIRTAG == pair_tags[1]:
            return f"{unit.fq2}"
    else:
        return f"{unit.fq1}"


def fastqc_trim_inputs(wildcards):
    if is_paired_end(wildcards.SAMPLE):
        if wildcards.PAIRTAG == pair_tags[0]:
            return "results/trim/fastq/{SAMPLE}_{UNIT}_R1.fastq.gz"
        elif wildcards.PAIRTAG == pair_tags[1]:
            return "results/trim/fastq/{SAMPLE}_{UNIT}_R2.fastq.gz"
    else:
        return "results/trim/fastq/{SAMPLE}_{UNIT}_R0.fastq.gz"


def trim_inputs(wildcards):
    unit = units.loc[wildcards.SAMPLE, wildcards.UNIT]
    if is_paired_end(wildcards.SAMPLE):
        return {"sample": [f"{unit.fq1}", f"{unit.fq2}"]}
    else:
        return {"sample": [f"{unit.fq1}"]}


def align_inputs(wildcards):
    if is_paired_end(wildcards.SAMPLE):
        return {
            "fq1": "results/trim/fastq/{SAMPLE}_{UNIT}_R1.fastq.gz",
            "fq2": "results/trim/fastq/{SAMPLE}_{UNIT}_R2.fastq.gz",
        }
    else:
        return {"fq1": "results/trim/fastq/{SAMPLE}_{UNIT}_R0.fastq.gz"}


def mark_duplicates_inputs(wildcards):
    sample_units = units.loc[wildcards.SAMPLE]
    if config["umi_dedup"]:
        if is_paired_end(wildcards.SAMPLE):
            return {
                "bams": expand(
                    "results/group_umis/bam/{{SAMPLE}}_{UNIT}_paired_end.bam",
                    UNIT=sample_units["unit"],
                ),
                "bams_idx": expand(
                    "results/group_umis/bam/{{SAMPLE}}_{UNIT}_paired_end.bam.bai",
                    UNIT=sample_units["unit"],
                ),
            }
        else:
            return {
                "bams": expand(
                    "results/group_umis/bam/{{SAMPLE}}_{UNIT}_single_end.bam",
                    UNIT=sample_units["unit"],
                ),
                "bams_idx": expand(
                    "results/group_umis/bam/{{SAMPLE}}_{UNIT}_single_end.bam.bai",
                    UNIT=sample_units["unit"],
                ),
            }
    else:
        return {
            "bams": expand(
                "results/assign_read_groups/bam/{{SAMPLE}}_{UNIT}.bam",
                UNIT=sample_units["unit"],
            ),
            "bams_idx": expand(
                "results/assign_read_groups/bam/{{SAMPLE}}_{UNIT}.bam.bai",
                UNIT=sample_units["unit"],
            ),
        }


def split_n_cigar_reads_inputs(wildcards):
    if config["umi_dedup"]:
        return {"bam": "results/mark_duplicates/bam/{SAMPLE}_umi.bam"}
    else:
        return {"bam": "results/mark_duplicates/bam/{SAMPLE}.bam"}


def trim_md5_inputs():
    inputs = []
    for sample in samples["sample"]:
        sample_units = units.loc[sample]
        inputs.extend(
            expand(
                "results/trim/fastq/{SAMPLE}_{UNIT}_{PAIR}.fastq.gz",
                SAMPLE=sample_units["sample"],
                UNIT=sample_units["unit"],
                PAIR=pair_tags,
            )
        )
    return inputs


def align_md5_inputs():
    inputs = []
    for sample in samples["sample"]:
        sample_units = units.loc[sample]
        inputs.extend(
            expand(
                "results/align/bam/{SAMPLE}_{UNIT}.bam",
                SAMPLE=sample_units["sample"],
                UNIT=sample_units["unit"],
            )
        )
    return inputs


def assign_read_groups_md5_inputs():
    inputs = []
    for sample in samples["sample"]:
        sample_units = units.loc[sample]
        inputs.extend(
            expand(
                "results/assign_read_groups/bam/{SAMPLE}_{UNIT}.bam",
                SAMPLE=sample_units["sample"],
                UNIT=sample_units["unit"],
            )
        )
    return inputs


def group_umis_md5_inputs():
    inputs = []
    for sample in samples["sample"]:
        sample_units = units.loc[sample]
        if is_paired_end(sample):
            inputs.extend(
                expand(
                    "results/group_umis/bam/{SAMPLE}_{UNIT}_paired_end.bam",
                    SAMPLE=sample_units["sample"],
                    UNIT=sample_units["unit"],
                )
            )
        else:
            inputs.extend(
                expand(
                    "results/group_umis/bam/{SAMPLE}_{UNIT}_single_end.bam",
                    SAMPLE=sample_units["sample"],
                    UNIT=sample_units["unit"],
                )
            )
    return inputs


def mark_duplicates_md5_inputs():
    if config["umi_dedup"]:
        return expand(
            "results/mark_duplicates/bam/{SAMPLE}_umi.bam",
            SAMPLE=samples["sample"],
        )
    else:
        return expand(
            "results/mark_duplicates/bam/{SAMPLE}.bam",
            SAMPLE=samples["sample"],
        )


####
## Workflow output files (Rule all inputs)
####


def workflow_outputs():
    """
    Returns all file endpoints for the workflow
    """

    outputs = []

    ## FastQC outputs
    for sample in samples["sample"]:
        sample_units = units.loc[sample]
        ## Raw
        outputs.extend(
            expand(
                "results/raw_data/FastQC/{SAMPLE}_{UNIT}_{PAIRTAG}_fastqc.{EXT}",
                SAMPLE=sample,
                UNIT=sample_units["unit"],
                PAIRTAG=pair_tags,
                EXT=["html", "zip"],
            )
        )
        ## Trim
        outputs.extend(
            expand(
                "results/trim/FastQC/{SAMPLE}_{UNIT}_{PAIRTAG}_fastqc.{EXT}",
                SAMPLE=sample,
                UNIT=sample_units["unit"],
                PAIRTAG=pair_tags,
                EXT=["html", "zip"],
            )
        )
        ## Align
        outputs.extend(
            expand(
                "results/align/FastQC/{SAMPLE}_{UNIT}_fastqc.{EXT}",
                SAMPLE=sample,
                UNIT=sample_units["unit"],
                EXT=["html", "zip"],
            )
        )

    ## md5sums
    outputs.extend(
        [
            "results/trim/fastq/md5.txt",
            "results/align/bam/md5.txt",
            "results/assign_read_groups/bam/md5.txt",
            "results/mark_duplicates/bam/md5.txt",
            "results/split_n_cigar_reads/bam/md5.txt",
            "results/bqsr/bam/md5.txt",
            "results/variation/gvcf/md5.txt",
            "results/variation/vcf/md5.txt",
        ]
    )
    if config["umi_dedup"]:
        outputs.append("results/group_umis/bam/md5.txt")

    ## Processed variants
    outputs.append("results/variation/vcf/snvs.final.vcf")

    return outputs
