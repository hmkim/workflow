library(dplyr)
library(readr)
args = commandArgs(trailingOnly = TRUE)
in_fn = args[1]
out_fn = args[2]

t = read_delim(in_fn, delim=" ")

x = t %>% filter(mapq > 20) %>%
  group_by(read_name, chrom, pos, orientation, insertion_end, seqcode) %>%
  distinct() %>%
  group_by(read_name) %>% mutate(read_count=n()) %>% filter(mapq == max(mapq)) %>%
  group_by(chrom, pos, orientation, insertion_end, seqcode) %>% summarise(count=n())
x = x[, c("chrom", "pos", "pos", "orientation", "insertion_end", "seqcode", "count")]

write.table(x, file=out_fn, row.names=FALSE, quote=FALSE, col.names=FALSE,
            sep="\t")
