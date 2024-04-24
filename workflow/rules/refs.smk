rule genome_get:
    output:
        temp("resources/genome.fa"),
    params:
        species=config["ref"]["species"],
        datatype="dna",
        build=config["ref"]["build"],
        release=config["ref"]["release"],
    wrapper:
        "v3.7.0/bio/reference/ensembl-sequence"

rule genome_index:
    input:
        "resources/genome.fa",
    output:
        temp("resources/genome.fa.fai"),
    params:
        extra="",
    wrapper:
        "v3.7.0/bio/samtools/faidx"

rule genome_dict:
    input:
        "resources/genome.fa",
    output:
        temp("resources/genome.dict"),
    params:
        extra="",
    wrapper:
        "v3.7.0/bio/picard/createsequencedictionary"

rule annotation_get:
    output:
        temp("resources/annotation.gtf"),
    params:
        species=config["ref"]["species"],
        build=config["ref"]["build"],
        release=config["ref"]["release"],
        flavor="",
    wrapper:
        "v3.7.0/bio/reference/ensembl-annotation"

rule star_index:
    input:
        fasta="resources/genome.fa",
        gtf="resources/annotation.gtf",
    output:
        temp(directory("resources/genome")),
    params:
        sjdbOverhang=int(config["read_length"]) - 1,
        extra="",
    wrapper:
        "v3.7.0/bio/star/index"

rule known_variants_get:
    input:
        fai="resources/genome.fa.fai"
    output:
        vcf=temp("resources/known_variants.vcf.gz"),
    params:
        species=config["ref"]["species"],
        build=config["ref"]["build"],
        release=config["ref"]["release"],
        type="all",
    wrapper:
        "v3.7.0/bio/reference/ensembl-variation"

rule known_variants_index:
    input:
        "resources/known_variants.vcf.gz",
    output:
        temp("resources/known_variants.vcf.gz.tbi"),
    params:
        "-p vcf",
    wrapper:
        "v3.7.0/bio/tabix/index"