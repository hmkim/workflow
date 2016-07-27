#!/usr/bin/env Rscript

install.packages(c("stringr", "data.table"))

# Load and sum all quant files in quant.files vectortors.
# Returns a single vector of TPM values with names as
# transcript names.
sum_abundances <- function(abundance_files) {
    abundances <- NULL
    for (i in 1:length(abundance_files)) {
        sample <- names(abundance_files)[i]
        a <- read.table(quant.files[i], sep="\t", header=T, stringsAsFactors=FALSE)
        if (i == 1) {
            abundances <- a$tpm
            names(abundances) <- a$target_id
        }
        else {
            abundances += a$tpm
        }
    }
    abundances
}

filter_transcripts <- function(abundance_files, cutoff=0.1, tx2gene=NULL) {
    abundances <- data.table::data.table(do.call(cbind, lapply(abundance_files, sum_abundances)))
    if (is.null(tx2gene)) {
        abundances$tx_id <- stringr::str_extract(rownames(abundances), 'ENSTR?[\\d\\.]+')
        abundances$gene_id <- stringr::str_extract(rownames(abundances), 'ENSGR?[\\d\\.]+')
    }
    else {
        if (typeof(tx2gene) == "character") {
            tx2gene <- read.table(tx2gene, sep="\t", header=TRUE, stringsAsFactors=FALSE)
        }
        colnames(tx2gene) <- c("tx_id", "gene_id")
        m <- match(rownames(abundances), tx2gene$tx_id)
        if (any(is.na(m))) {
            stop("One or more transcript IDs missing from tx2gene mapping")
        }
        abundances <- cbind(abundances, data.table::as.data.table(tx2gene[m,]))
    }
    abundances[, rel_abundance := tpm / sum(tpm), gene_id]
    abundances[rel_abundance < cutoff, tx_id]
}

quant_files <- read.table($abundanceFile, sep="\t", header=F, stringsAsFactors=F)
quant_files <- tapply(quant_files[,2], quant_files[,1], identity, simplify=FALSE)
exclude_tx_ids <- filter_transcripts(quant_files, 
    as.numeric(${params.cutoff}), ${params.tx2gene_file})
writeLines(exclude_tx_ids, $excludeFile)
