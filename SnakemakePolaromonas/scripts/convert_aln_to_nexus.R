library(seqinr)
library(ape)
data = read.fasta(snakemake@input[['aln']])
names(data) <- gsub("-", "_", names(data))
write.nexus.data(
  data,
  file = snakemake@output[['nex']],
  format = "dna",
  interleaved = FALSE
)
