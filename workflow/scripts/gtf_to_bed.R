gtf <- read.table(snakemake@input[["gtf"]], sep="\t")
gtf <- subset(gtf, V3 == "exon")
df <- data.frame(chrom=gtf[,'V1'], start=gtf[,'V4'], end=gtf[,'V5'])
## GTF is 1-indexed but BED is 0-indexed
## Additionally, the 'end' coordinate for BED is not included in the feature
## e.g. GTF 1-100 == BED 0-100
df[["start"]] <- df[["start"]] - 1
write.table(
    df,
    snakemake@output[["bed"]],
    quote = FALSE,
    sep="\t",
    col.names = FALSE,
    row.names = FALSE
)
