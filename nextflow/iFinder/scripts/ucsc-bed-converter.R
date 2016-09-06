x = "/Users/rory/cache/li_hiv/metadata/hg19_gene_features.bed"
bed = read.table(x, header=FALSE, sep="\t", stringsAsFactors=FALSE)

get_transcript_id = function(id) {
    transcript_id = strsplit(id, "_")[[1]][1]
    return(transcript_id)
}
    
get_feature = function(id) {
    feature = strsplit(id, "_")[[1]][2]
    return(feature)
}

bed$feature = unlist(lapply(bed$V4, get_feature))
bed$transcript_id = unlist(lapply(bed$V4, get_transcript_id))

bed = bed[, c("V1", "V2", "V3", "V6", "transcript_id", "feature")]

library(biomaRt)
human = useMart("ensembl", dataset="hsapiens_gene_ensembl")
conversions = getBM(attributes=c("ensembl_gene_id", "ensembl_transcript_id",
                        "hgnc_symbol"), mart=human)

library(dplyr)

collapse = function(x) {
    l = list(x)
    return(paste(l[[1]], collapse=","))
}
collapsed = conversions %>% group_by(ensembl_transcript_id) %>%
    summarise(conv=collapse(hgnc_symbol), id=collapse(ensembl_gene_id))
rownames(collapsed) = collapsed$ensembl_transcript_id

m = merge(bed, collapsed, by.x="transcript_id", by.y="ensembl_transcript_id", all.x=TRUE)
#bed$symbol = collapsed[bed$transcript_id,]$conv

m = m[, c("V1", "V2", "V3", "conv", "feature", "id", "V6")]
write.table(m, file="features.bed", col.names=FALSE, sep="\t", quote=FALSE, row.names=FALSE)
