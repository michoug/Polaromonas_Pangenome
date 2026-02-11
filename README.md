# Polaromonas Pangenome Analysis

[![DOI](https://zenodo.org/badge/1155438049.svg)](https://doi.org/10.5281/zenodo.18611200)

This repository contains a comprehensive pangenome analysis pipeline for *Polaromonas* genomes using R and the `targets` package for reproducible workflow management.

## Overview

This project performs statistical analysis and visualization of *Polaromonas* pangenome data, including:
- Genome metadata analysis and quality assessment
- Pangenome structure analysis (core, shell, and cloud genes)
- Phylogenetic tree visualization and ancestral state reconstruction
- Gene family evolution analysis using duplication-transfer-loss-origination (DTLO) events
- Microbial trait prediction and enrichment analysis
- Geographic distribution mapping
- Comparative genomics and NMDS ordination

## Project Structure

### R Scripts (`R/` folder)

The `R/` directory contains modular R scripts that define all custom functions used in the analysis pipeline:

- **`packages.R`**: Loads all required R packages for the analysis (tidyverse, ggplot2, phylogenetic tools, etc.)
- **`functions.R`**: Core data processing functions for cleaning, transforming, and analyzing pangenome data
- **`plots.R`**: Visualization functions for generating figures and plots
- **`diverse.R`**: Additional utility functions (currently commented out in the pipeline)

### Workflow Pipeline (`_targets.R`)

The `_targets.R` file defines the entire analysis workflow using the [`targets`](https://docs.ropensci.org/targets/) package, which provides:

- **Reproducible pipeline**: Automatically tracks dependencies between analysis steps
- **Efficient computation**: Only re-runs outdated targets when data or code changes
- **Parallel processing**: Utilizes multiple CPU cores for faster execution
- **Organized outputs**: Generates figures, tables, and intermediate data objects

The pipeline includes:
1. **Data loading**: Reads genome statistics, quality metrics, pangenome data, phylogenetic trees, and functional annotations
2. **Data processing**: Cleans and combines metadata, performs statistical analyses
3. **Ordination analysis**: NMDS ordination on persistent, shell, and cloud genes
4. **Phylogenetic analysis**: Processes phylogenetic trees and DTL events
5. **Functional analysis**: Analyzes KEGG pathways, modules, and microbial traits
6. **Visualization**: Generates publication-ready figures (maps, heatmaps, trees, ordinations)

## Requirements

### R Environment

This project uses `renv` for R package management to ensure reproducibility:

```bash
# Restore the R environment
R -e "renv::restore()"
```

### Key R Packages

- **Workflow**: [targets](https://docs.ropensci.org/targets/), [tarchetypes](https://docs.ropensci.org/tarchetypes/), [crew](https://docs.ropensci.org/crew/)
- **Data manipulation**: [tidyverse](https://www.tidyverse.org/), [vroom](https://vroom.r-lib.org/), [janitor](https://sfirke.github.io/janitor/)
- **Phylogenetics**: [ape](https://cran.r-project.org/web/packages/ape/index.html), [ggtree](https://bioconductor.org/packages/release/bioc/html/ggtree.html), [treeio](https://bioconductor.org/packages/release/bioc/html/treeio.html), [tidytree](https://bioconductor.org/packages/release/bioc/html/tidytree.html)
- **Visualization**: [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html), [patchwork](https://cran.r-project.org/web/packages/patchwork/index.html), [gt](https://cran.r-project.org/web/packages/gt/index.html), [gtExtras](https://cran.r-project.org/web/packages/gtExtras/index.html)
- **Statistics**: [vegan](https://cran.r-project.org/web/packages/vegan/index.html) (for NMDS and ANOSIM)
- **Geospatial**: [sf](https://cran.r-project.org/web/packages/sf/index.html), [rnaturalearth](https://cran.r-project.org/web/packages/rnaturalearth/index.html)

## Usage

### Running the Analysis Pipeline

Execute the complete workflow using the `targets` package:

```r
# Load targets library
library(targets)
library(tarchetypes)
library(crew)
# Run with parallel processing (configured in _targets.R)
tar_make(reporter = "balanced")
```

### Viewing Pipeline Results

```r
# List all targets
tar_manifest()

# Load specific results
tar_load(figure_phylogenetic_tree)
tar_load(genome_metadata_comb)

# Read specific targets
tar_read(figure_map_world)
```

### Generating Specific Outputs

To regenerate specific figures or analyses:

```r
# Regenerate a single target
tar_make(figure_phylogenetic_tree)

# Regenerate all figures
tar_make(starts_with("figure_"))
```

## Output Structure

The pipeline generates outputs in the `Figures/` directory:

- **Figure 1**: Geographic distribution maps and genome comparisons
- **Figure 2**: NMDS ordination of persistent genes
- **Figure 3**: Phylogenetic tree with ancestral state reconstruction
- **Figure 4**: Heatmap of DTL events and functional traits
- **Figure S1-S6**: Supplementary figures (rarefaction curves, prevalence plots, etc.)
- **Table S1-S2**: Metadata and genome statistics tables

## Data Requirements

All the data that the pipeline expects are found in the `rawData/` directory:

- Genome statistics and quality metrics (CheckM2 results)
- Pangenome data (PPanGGOLiN output)
- Gene presence/absence matrices
- Phylogenetic trees
- Functional annotations (eggNOG, KEGG)
- Microbial trait predictions
- Geographic metadata

## Upstream Workflow

For genome preprocessing and pangenome construction, see the **`SnakemakePolaromonas/`** folder. This folder contains a Snakemake pipeline for genome annotation (Bakta), pangenome analysis (PPanGGOLiN), and phylogenetic tree construction. Refer to its README for details.

## Citation

If you use this pipeline or data, please cite the associated publication (details to be added).

## License

MIT License.

## Contact

Grégoire Michoud
