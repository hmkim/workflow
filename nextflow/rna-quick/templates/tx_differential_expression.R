source("http://bioconductor.org/biocLite.R")
biocLite("devtools")
biocLite("pachterlab/sleuth")

load("$sleuthObj")

name <- "$modelName"
formula <- $modelFormula
betas <- "$modelBetas"
betas <- unlist(strsplit(betas, ","))

sleuth_data <- sleuth::sleuth_fit(sleuth_data, formula=formula, fit_name=name)

for (i in 1:length(betas)) {
    b <- strsplit(betas[[i]], "=")[[1]]
    stopifnot(b[1] %in% colnames(sleuth_data$obs_raw))
    if (length(b) == 1) {
        val <- sleuth_data$obs_raw[,b[1]][1]
    }
    else {
        val <- b[2]
    }
    beta <- paste0(beta[[1]], val)
    sleuth_data <- sleuth::sleuth_wt(sleuth_data, which_beta=beta, which_mod=name)
}

save(sleuth_data, file="$resultFile")
