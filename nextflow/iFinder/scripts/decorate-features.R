library(dplyr)
library(tools)
args = commandArgs(trailingOnly=TRUE)
gene_fn = args[1]
virus_fn = args[2]
out_fn = paste(dirname(gene_fn), "/", strsplit(basename(gene_fn), "-gene-features")[[1]][1],
    "-decorated.bed", sep="")
gene_features = read.table(gene_fn, header=FALSE, sep="\t", stringsAsFactors=FALSE)
virus_features = read.table(virus_fn, header=FALSE, sep="\t", stringsAsFactors=FALSE)

colnames(gene_features) = c("chrom", "start", "end", "sample", "virus_start", "orientation",
            "count", "gene", "feature", "id", "strand")
colnames(virus_features) = c("virus_chrom", "virus_start", "virus_end", "virus_feature")
virus_features = virus_features %>% group_by(virus_chrom, virus_start, virus_end) %>%
    summarise_each(funs(toString(unique(.))))
m = gene_features %>% left_join(virus_features, by="virus_start")
m = m[complete.cases(m),]
remap_3ltr = m$virus_start >= 9085 & m$virus_start <= 9540
m$virus_feature[remap_3ltr] = "5-ltr"
m$virus_start[remap_3ltr] = m[remap_3ltr,]$virus_start - 9085
m$virus_end[remap_3ltr] = m[remap_3ltr,]$virus_start + 1

remap_ru5 = m$virus_start > 456 & m$virus_start < 634
m$virus_feature[remap_ru5] = "nef, 3-ltr"
m$virus_start[remap_ru5] = m[remap_ru5,]$virus_start + 9085
m$virus_end[remap_ru5] = m[remap_ru5,]$virus_start + 1
m = m %>% group_by(chrom, start, end, sample, virus_start, orientation,  gene, id,
    feature, strand, virus_chrom, virus_end, virus_feature) %>% summarise(count=sum(count))
m = m %>% group_by(chrom, start, end, sample, virus_start, orientation, virus_chrom,
    virus_end, virus_feature) %>% summarise_each(funs(toString(unique(.))))
write.table(m, file=out_fn, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")


