library(ape)
library(tidyverse)


print(snakemake@input[['nex']])

run1 = read.nexus(snakemake@input[['nex']])
run2 = read.nexus(gsub(".nex.run1.t", ".nex.run2.t", snakemake@input[['nex']]))
print("writing")
write.tree(
  c(run1[502:1001], run2[502:1001]),
  file = snakemake@output[['newick']]
)

mapping <- data.frame((run1[1000]$tip.label), (run1[1000]$tip.label))
mapping$new <- paste0(mapping$tip.label, sep = ":", mapping$tip.label.1)

write.table(
  mapping[, "new"],
  file = snakemake@output[['txt']],
  col.names = F,
  row.names = F,
  quote = F
)
