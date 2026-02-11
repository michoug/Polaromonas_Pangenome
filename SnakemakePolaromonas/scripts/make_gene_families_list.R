library(ape)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

file_list <- list.files(path = "args[1]")

names <- gsub(".newick", "", file_list)

print(dim(names))

families_list <- c("[FAMILIES]")

for (i in names) {
  tmp <- c(
    paste0("- ", i),
    paste0("gene_tree = gene_trees/", i, ".newick"),
    paste0("mapping = mapping_files/", i, ".txt")
  )
  families_list <- c(families_list, tmp)
}

write(families_list, file = "args[2]")
