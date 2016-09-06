library(dplyr)
library(readr)
args = commandArgs(trailingOnly = TRUE)
in_fn = args[1]
out_fn = args[2]

t = read_delim(in_fn, delim=" ")

left = subset(t, orientation == "5prime")
right = subset(t, orientation == "3prime")

left = left %>% filter(mapq > 20) %>%
  group_by(read_name, chrom, pos, orientation, insertion_end, seqcode, virus_pos,
           mapq) %>%
  distinct() %>%
  group_by(read_name) %>% mutate(read_count=n()) %>% filter(mapq == max(mapq)) %>%
  mutate(nmult=n()) %>%
  group_by(chrom, pos, orientation, insertion_end, seqcode, virus_pos) %>%
  summarise(count=n(), amult=mean(nmult))
left = left[, c("chrom", "pos", "pos", "orientation", "insertion_end",
                "seqcode", "count", "virus_pos", "amult")]
colnames(left) = c("chrom", "start", "end", "orientation", "insertion_end", "seqcode",
                   "count", "virus_pos", "amult")

right = right %>% filter(mapq > 20) %>%
  group_by(read_name, chrom, pos, orientation, insertion_end, seqcode, virus_end,
           mapq) %>%
  distinct() %>%
  group_by(read_name) %>% mutate(read_count=n()) %>% filter(mapq == max(mapq)) %>%
  mutate(nmult=n()) %>%
  group_by(chrom, pos, orientation, insertion_end, seqcode, virus_end) %>%
  summarise(count=n(), amult=mean(nmult))
right = right[, c("chrom", "pos", "pos", "orientation", "insertion_end",
                  "seqcode", "count", "virus_end", "amult")]
colnames(right) = c("chrom", "start", "end", "orientation", "insertion_end",
                    "seqcode", "count", "virus_pos", "amult")
if(nrow(right) > 0) {
  all = rbind(left, right)
} else {
  all = left
}
write.table(all, file=out_fn, row.names=FALSE, quote=FALSE, col.names=FALSE,
            sep="\t")
