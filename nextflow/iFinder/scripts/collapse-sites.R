library(dplyr)
library(tools)
args = commandArgs(trailingOnly = TRUE)
fn = args[1]
outfn = paste(file_path_sans_ext(fn), "-sites.bed", sep="")
t = read.table(fn, header=FALSE, sep= "\t", stringsAsFactors=FALSE)
colnames(t) = c("chrom", "start", "end", "rid", "sample", "virus_start", "orientation")
grouped = t %>% group_by(chrom, start, end, rid, sample, virus_start, orientation) %>%
    summarise_each(funs(toString(Filter(function(x) x != "", unique(.)))))

sites = grouped %>% group_by(chrom, start, end, sample, virus_start,
                             orientation) %>% summarise(count=n())
sites$start = sites$start + 1
sites$end = sites$end + 1
write.table(sites, file=outfn, col.names=FALSE, row.names=FALSE,
            quote=FALSE, sep="\t")
