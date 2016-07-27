source("http://bioconductor.org/biocLite.R")
biocLite("DESeq2")
biocLite("tximport")

load("$sleuthObj")

name <- "$modelName"
formula <- $modelFormula
betas <- "$modelBetas"
betas <- unlist(strsplit(betas, ","))
threads <- $cpus

get_eff_len <- function(sleuth_data) {
    s_data <- as.data.frame(sleuth_data[["obs_norm"]]) %>%
        select_("target_id", "sample", "eff_len") %>%
        tidyr::spread_("sample", which_units)
    rownames(s_data) <- s_data$target_id
    s_data$target_id <- NULL
    s_data <- as.matrix(s_data)
    s_data[, sample(1:ncol(s_data))]
}

txi <- list(
    abundance=sleuth_to_matrix(sleuth_data, "obs_raw", "tpm")$data,
    counts=sleuth_to_matrix(sleuth_data, "obs_raw", "est_counts")$data
    length=get_eff_len(sleuth_data),
    countsFromAbundance=NULL)
)
gene_data <- tximportData::summarizeToGene(txi, tx2gene, ignoreTxVersion, "no")

dds <- DESeq2::DESeqDataSetFromMatrix(
    countData=gene_data$counts, 
    colData=sleuth_data$sample_to_covariates,
    design=formula
)

if (threads == 1) {
    dds <- DESeq2::DESeq(dds)
}
else {
    biocLite("BiocParallel")
    dds <- DESeq2::DESeq(dds, parallel=TRUE, BPPARAM=BiocParallel::MulticoreParam(workers=threads))
}

deseq2_data <- results(dds)
save(deseq2_data, file="$resultFile")
