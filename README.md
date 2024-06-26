# Snakemake workflow for variant calling from RNA-seq data

This Snakemake workflow implements the [GATK best-practices workflow for RNA-seq short variant discovery](https://gatk.broadinstitute.org/hc/en-us/articles/360035531192-RNAseq-short-variant-discovery-SNPs-Indels).

Complete documentation coming when I finish my PhD... but the workflow is fully functional as far as I am aware.

It has been tested for compatibility with:
- Single end reads
- Paired end reads
- Multiple sequencing units for one or more samples
- Deduplication with or without UMIs

Configuration documentation is however available in `config/README.md`, and as comments in `config/config.yaml`.
