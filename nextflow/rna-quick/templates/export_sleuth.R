source("http://bioconductor.org/biocLite.R")
biocLite("devtools")
biocLite("pachterlab/sleuth")

libraries_file <- "$libraries"
abundances_file <- "$abundances"
full_model <- $fullModel
tx2gene_file <- "$tx2gene"

metadata <- read.table(libraries_file, sep="\t", header=T, stringsAsFactors=F)
abundance_files <- read.table(abundances_file, sep="\t", header=F, stringsAsFactors=F)
m <- match(metadata[,"UniqueID"], abundance_files[,1])
if (any(is.na(m))) {
    stop("One or more samples missing from list of abundance files")
}
metadata[,"path"] <- abundance_files[m, 2]

sleuth_data <- sleuth::sleuth_prep(metadata, as.formula(full_model))

if (file.exists(tx2gene_file)) {
    tx2gene <- read.table(tx2gene, sep="\t", header=T, stringsAsFactors=F)
    m <- match(rownames(sleuth_data[["obs_raw"]]), tx2gene[,"target_id"])
    if (any(is.na(m))) {
        stop(paste(tx2gene_file, "missing one or more transcript IDs"))
    }
}
else {
    install.packages(stringr)
    tx2gene <- data.frame(
        tx_id=stringr::str_extract(rownames(abundance_matrix), 'ENSTR?[\\d\\.]+'),
        gene_id=stringr::str_extract(rownames(abundance_matrix), 'ENSGR?[\\d\\.]+')
    )
}

sleuth_data[["target_mapping"]] <- tx2gene

save(sleuth_data, file="$sleuthObj")