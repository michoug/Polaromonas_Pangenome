library <- function(...) {
  packages <- as.character(match.call(expand.dots = FALSE)[[2]])
  suppressWarnings(suppressMessages(lapply(
    packages,
    base::library,
    character.only = TRUE
  )))
  return(invisible())
}

## load the packages

library(tidyverse)

library(ape)
library(fs)
library(ggimage)
library(ggnewscale)
library(ggpp)
library(ggrepel)
library(ggtree)
library(ggtreeExtra)
library(glue)
library(gt)
library(gtExtras)
library(janitor)
library(MonoPhy)
library(patchwork)
library(RColorBrewer)
library(rnaturalearth)
library(rnaturalearthdata)
library(scales)
library(sf)
library(tidytree)
library(treeio)
library(vegan)
library(vroom)
library(webshot2)
