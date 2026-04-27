theme_all <- ggplot2::theme(
  text = ggplot2::element_text(size = 12, colour = "black"),
  strip.text = ggplot2::element_text(size = 14, colour = "black"),
  strip.text.x.top = ggplot2::element_text(angle = 90),
  axis.text = ggplot2::element_text(size = 12, colour = "black"),
  strip.text.y.left = ggplot2::element_text(angle = 0),
  strip.placement = "outside",
  panel.background = ggplot2::element_blank(),
  legend.title = ggplot2::element_text(
    size = 12,
    face = "bold",
    colour = "black"
  ),
  legend.text = ggplot2::element_text(size = 12, colour = "black"),
  legend.position = 'top',
  legend.box = 'vertical',
  legend.key = ggplot2::element_blank(),
  legend.margin = ggplot2::margin(b = -5),
)

# Creates a small pie chart image saved to disk.
# Used to overlay environment composition on map figures.
create_pie_image <- function(data, filename, colors_samples) {
  pie <- ggplot(data, aes(x = "", y = value, fill = env_label_good)) +
    geom_bar(stat = "identity", width = 1, show.legend = FALSE) +
    coord_polar("y") +
    theme_void() +
    scale_fill_manual(values = colors_samples)

  ggsave(
    filename,
    plot = pie,
    width = 2,
    height = 2,
    dpi = 100,
    bg = "transparent",
    create.dir = TRUE
  )
}

# Plots a world map of Polaromonas relative abundances from MicrobeAtlas.
# Splits Europe into an inset panel for readability.
plot_map_atlas <- function(atlas_map) {
  world <- ne_countries(scale = "medium", returnclass = "sf")

  colors_atlas <- c(brewer.pal(12, "Paired"), "black")
  colors_atlas <- gsub("#FFFF99", "lightgoldenrod3", colors_atlas)
  names(colors_atlas) <- levels(atlas_map$env_good)

  atlas_map_no_eur <- atlas_map |>
    filter(
      !(parsedLon >= -10.3 &
        parsedLon <= 26.5 &
        parsedLat >= 35.6 &
        parsedLat <= 71.4)
    )

  atlas_map_eur <- atlas_map |>
    filter(
      (parsedLon >= -10.3 &
        parsedLon <= 26.5 &
        parsedLat >= 35.6 &
        parsedLat <= 71.4)
    )

  p1 <- ggplot() +
    geom_sf(data = world, fill = "white") +
    geom_point(
      data = atlas_map_no_eur,
      aes(
        x = parsedLon,
        y = parsedLat,
        color = env_good,
        size = percAbund
      )
    ) +
    scale_color_manual(values = colors_atlas, guide = "none") +
    scale_size(limits = c(0, 12), breaks = c(0.5, 1, 2, 5, 10)) +
    coord_sf(
      xlim = c(-160, 170),
      ylim = c(-80, 85),
      expand = FALSE
    ) +
    geom_rect(
      data = world,
      xmin = -10.3,
      xmax = 26.5,
      ymin = 35.6,
      ymax = 71.4,
      fill = NA,
      colour = "black",
      linewidth = 1
    ) +
    labs(
      color = "Environment",
      x = NULL,
      y = NULL,
      size = "Relative\nAbundance\n(%)"
    ) +
    theme_light() +
    theme(
      axis.text = element_text(color = "black"),
      text = element_text(size = 12, colour = "black"),
      legend.text = element_text(size = 12, colour = "black")
    )

  p2 <- ggplot() +
    geom_sf(data = world, fill = "white") +
    geom_point(
      data = atlas_map_eur,
      aes(
        x = parsedLon,
        y = parsedLat,
        color = env_good,
        size = percAbund
      )
    ) +
    scale_color_manual(values = colors_atlas) +
    scale_size(limits = c(0, 12), breaks = c(0.5, 1, 2, 5, 10)) +
    coord_sf(
      xlim = c(-10.3, 26.5),
      ylim = c(35.6, 71.4),
      expand = FALSE
    ) +
    labs(
      color = "Environment",
      x = NULL,
      y = NULL,
      size = "Relative\nAbundance\n(%)"
    ) +
    theme_light() +
    theme(
      axis.text = element_text(color = "black"),
      text = element_text(size = 12, colour = "black"),
      legend.text = element_text(size = 12, colour = "black")
    )

  p3 <- wrap_plots(p1, p2, nrow = 1, guides = "collect")
  p3
}

# Prepares spatial data (centroids, pie images, genome counts) for map plotting.
# Generates per-country pie chart PNGs and merges with shapefile centroids.
prepare_map_data <- function(
  genome_metadata_clean,
  colors_samples,
  region = c("world", "europe")
) {
  region <- match.arg(region)

  # Count genomes per country
  numbgenome <- genome_metadata_clean |>
    group_by(country) |>
    summarise(numb = n(), .groups = "drop")

  # Load world shapefile
  world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

  # Filter countries present in your dataset
  world_select <- world |>
    filter(name %in% genome_metadata_clean$country) |>
    left_join(numbgenome, by = c("name" = "country")) |>
    mutate(name = factor(name, levels = sort(name)))

  # Pie data
  pie_data <- genome_metadata_clean |>
    group_by(country, env_label_good) |>
    summarise(value = n(), .groups = "drop")

  # Centroids
  centroids <- st_point_on_surface(world_select)
  centroids_coords <- as.data.frame(st_coordinates(centroids))
  centroids_coords$name <- centroids$name

  # Selected countries
  selected_countries <- sort(world_select$name)

  # Create a proper mapping: country -> pie image file
  image_df <- data.frame(
    name = selected_countries,
    image = paste0(
      "Figures/temp/pie_",
      gsub(" ", "_", selected_countries),
      ".png"
    ),
    stringsAsFactors = FALSE
  )

  # Create pie images
  for (i in seq_along(selected_countries)) {
    country_data <- pie_data |> filter(country == selected_countries[i])
    create_pie_image(country_data, image_df$image[i], colors_samples)
  }

  # Merge centroids with images and genome counts
  centroids_df <- centroids_coords |>
    left_join(image_df, by = "name") |> # Correct mapping of image -> country
    left_join(numbgenome, by = c("name" = "country"))

  # Filter for Europe if needed
  european_countries <- c(
    "Austria",
    "Czechia",
    "France",
    "Germany",
    "Switzerland",
    "Norway",
    "Sweden",
    "Spain",
    "Italy",
    "United Kingdom"
  )

  if (region == "europe") {
    centroids_df <- centroids_df |> filter(name %in% european_countries)
  } else if (region == "world") {
    centroids_df <- centroids_df |> filter(!name %in% european_countries)
  }

  # Return list for plotting
  list(
    world = world,
    world_select = world_select,
    centroids = centroids_df,
    pie_data = pie_data
  )
}

# Renders the genome distribution map with pie overlays and country labels.
# Adjusts coordinate extent based on the requested region (world or Europe).
plot_map <- function(map_data, colors_samples, region = c("world", "europe")) {
  region <- match.arg(region)

  # Fake pie for legend (all environments)
  pie_fake <- data.frame(
    name = factor(names(colors_samples), levels = names(colors_samples)),
    X = 1,
    Y = 1
  )

  centroids_df <- map_data$centroids |>
    mutate(
      nudge_X = case_when(
        name == "Norway" ~ 8.5 - X,
        name == "New Zealand" ~ 170 - X,
        name == "Austria" ~ 4,
        TRUE ~ 0
      ),
      nudge_Y = case_when(
        name == "Nepal" ~ -2,
        name == "Norway" ~ -19,
        TRUE ~ 0
      ),
      X_nudged = X + nudge_X,
      Y_nudged = Y + nudge_Y
    )

  # Base map
  p <- ggplot() +
    geom_sf(data = map_data$world, fill = "white") +
    geom_sf(data = map_data$world_select, aes(fill = numb), color = "black") +
    scale_fill_gradient2() +
    geom_image(
      data = centroids_df,
      aes(x = X_nudged, y = Y_nudged, image = image),
      size = 0.075
    ) +
    geom_label_repel(
      data = centroids_df,
      aes(label = name, x = X_nudged, y = Y_nudged),
      box.padding = 0.5,
      point.padding = 0.3,
      nudge_x = case_when(centroids_df$name == "Austria" ~ 4, .default = 2),
      nudge_y = case_when(
        centroids_df$name == "Italy" ~ -2,
        .default = 2
      ),
      size = 3,
      fontface = if_else(centroids_df$numb > 5, "bold", "plain")
    ) +
    labs(fill = "Number of\nGenomes", x = NULL, y = NULL) +
    new_scale_fill() +
    geom_col(data = pie_fake, aes(x = 0, y = 0, fill = name)) +
    scale_fill_manual(values = colors_samples) +
    labs(fill = if (region == "world") "Environment" else "Environment") +
    theme_all +
    theme_light() +
    theme(axis.text = element_text(color = "black"))

  # Optionally adjust coordinates
  if (region == "world") {
    p <- p + coord_sf(xlim = c(-130, 180), ylim = c(-90, 90), expand = FALSE)
  } else {
    p <- p +
      coord_sf(xlim = c(-10.3, 26.5), ylim = c(35.6, 71.4), expand = TRUE)
  }

  p
}

# Builds a gt summary table of genome statistics per environment.
# Reports mean ± SD for length, GC, completeness, and contamination.
draw_table <- function(genome_statToPlot) {
  tab_all <- genome_statToPlot |>
    select(
      env_label_good,
      n,
      relative_length,
      GC,
      Completeness,
      Contamination
    ) |>
    mutate(relative_length = relative_length / 1e6) |>
    summarise(
      across(
        relative_length:Contamination,
        .fns = list(
          mean = function(x) mean(x, na.rm = TRUE),
          sd = function(x) sd(x, na.rm = TRUE)
        ),
        .names = "{.col}_{.fn}"
      )
    ) |>
    mutate(n = 283) |>
    mutate(env_label_good = "All")

  tab <- genome_statToPlot |>
    group_by(env_label_good) |>
    mutate(n = n()) |>
    ungroup() |>
    select(
      env_label_good,
      n,
      relative_length,
      GC,
      Completeness,
      Contamination
    ) |>
    mutate(relative_length = relative_length / 1e6) |>
    group_by(env_label_good, n) |>
    summarise(across(
      where(is.numeric),
      .fns = list(
        mean = function(x) mean(x, na.rm = TRUE),
        sd = function(x) sd(x, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    )) |>
    ungroup() |>
    as.data.frame() |>
    rbind(tab_all) |>
    gt() |>
    # tab_header(title = "Polaromonas Statistics", ) |>
    fmt_number(starts_with("relative")) |>
    fmt_number(
      starts_with("GC") |
        starts_with("Completeness") |
        starts_with("Contamination"),
      decimals = 1
    ) |>
    cols_merge_uncert(
      col_val = "relative_length_mean",
      col_uncert = "relative_length_sd"
    ) |>
    cols_merge_uncert(
      col_val = "Completeness_mean",
      col_uncert = "Completeness_sd"
    ) |>
    cols_merge_uncert(
      col_val = "Contamination_mean",
      col_uncert = "Contamination_sd"
    ) |>
    cols_merge_uncert(col_val = "GC_mean", col_uncert = "GC_sd") |>
    cols_label(
      env_label_good = "Environment",
      n = "Number of\ngenomes",
      relative_length_mean = "Relative length\n(Mbp)",
      GC_mean = "GC (%)",
      Completeness_mean = "Completeness\n(%)",
      Contamination_mean = "Contamination\n(%)"
    ) |>
    cols_align(align = "center", columns = everything())
  tab
}

# Plots violin + boxplot distributions of genome parameters per environment.
# Currently displays normalized genome length and GC content.
comparing_genomes_params <- function(genome_statToPlot, colors_samples) {
  plot_labels <- c(
    relative_length = "Normalized\nlength (Mbp)",
    GC = "GC content (%)"
  )

  p <- genome_statToPlot |>
    select("relative_length", "GC", env_label_good) |>
    mutate(relative_length = relative_length / 1e6) |>
    pivot_longer(cols = !env_label_good) |>
    mutate(
      name = factor(
        name,
        levels = (c(
          "relative_length",
          "GC"
        ))
      )
    ) |>
    ggplot(aes(x = env_label_good, y = value, fill = env_label_good)) +
    geom_violin() +
    geom_boxplot(width = 0.1, fill = "white") +
    facet_wrap(~name, scales = "free_y", labeller = as_labeller(plot_labels)) +
    theme_classic() +
    theme_all +
    theme(
      strip.text.x.top = element_text(angle = 0, size = 12),
      axis.text.x = element_text(
        angle = 45,
        vjust = 1,
        hjust = 1
      ),
      axis.title = element_blank()
    ) +
    labs(y = "", x = "") +
    scale_fill_manual(values = colors_samples, guide = "none")

  p
}

# Plots a 2D NMDS ordination with spider segments to environment centroids.
# Annotates the stress value on the figure.
plot_nmds <- function(nmds, map, colors_samples) {
  datfort <- as.data.frame(scores(nmds)$sites)

  datfort_sites <- datfort |>
    rownames_to_column(var = "label")

  nmds$stress <- ifelse(
    nmds$stress < 1e-2,
    scientific(nmds$stress, digits = 2),
    round(nmds$stress, digits = 2)
  )

  stress_value <- paste("Stress = ", nmds$stress, sep = "")

  yvalue = max(datfort_sites$NMDS2) -
    0.05 *
      (max(datfort_sites$NMDS2) - min(datfort_sites$NMDS2))
  xvalue = min(datfort_sites$NMDS1) +
    0.1 *
      (max(datfort_sites$NMDS1) - min(datfort_sites$NMDS1))

  NMDS_species1 <- nmds$species[, 1]
  NMDS_species2 <- nmds$species[, 2]
  NMDS_species <- data.frame(NMDS_species1, NMDS_species2) |>
    rownames_to_column("Modules")
  NMDS1 = nmds$points[, 1]
  NMDS2 = nmds$points[, 2]

  NMDS_data <- data.frame(NMDS1, NMDS2) |>
    rownames_to_column("Genome") |>
    mutate(Genome = gsub("-", "_", Genome)) |>
    left_join(map, join_by(Genome)) |>
    filter(env_label_good != "Other")

  cent <- aggregate(
    cbind(NMDS1, NMDS2) ~ env_label_good,
    data = NMDS_data,
    FUN = mean
  )

  segs <- merge(
    NMDS_data,
    setNames(cent, c('env_label_good', 'oNMDS1', 'oNMDS2')),
    by = 'env_label_good',
    sort = FALSE
  )

  p <- ggplot() +
    geom_point(
      data = NMDS_data,
      mapping = aes(x = NMDS1, y = NMDS2, colour = env_label_good),
      size = 2
    ) +
    geom_segment(
      data = segs,
      mapping = aes(
        x = NMDS1,
        y = NMDS2,
        colour = env_label_good,
        xend = oNMDS1,
        yend = oNMDS2
      )
    ) + # spiders
    geom_point(
      data = cent,
      mapping = aes(x = NMDS1, y = NMDS2, fill = env_label_good),
      colour = "black",
      shape = 21,
      size = 5,
      show.legend = FALSE
    ) +
    scale_fill_manual(values = colors_samples) +
    scale_color_manual(values = colors_samples) +
    annotate(
      "text",
      x = xvalue,
      y = yvalue,
      label = stress_value,
      size = 5
    ) +
    coord_equal() +
    labs(color = "Environment") +
    guides(fill = guide_legend(ncol = 1, byrow = TRUE)) +
    theme_classic() +
    theme_all

  p
}

# Draws a heatmap of enriched microtrait percentages across environments.
# Clusters traits by Euclidean distance and facets by trait category.
plot_heatmap_microtrait <- function(
  microtrait_to_plot,
  microtrait_enrich_sign
) {
  plot_data <- microtrait_to_plot |>
    filter(name %in% microtrait_enrich_sign$name) |>
    summarise(
      .by = c(env_label_good, name),
      per_genomes = sum(value) / n_distinct(Genome) * 100
    ) |>
    separate(
      name,
      into = c("trait_category1", "trait_category2", "trait_category3"),
      sep = "__",
      extra = "merge",
      remove = FALSE
    ) |>
    filter(env_label_good != "Other") |>
    mutate(
      trait_category3 = gsub(
        "(.*?):.*:(.*?)",
        "\\1 - \\2",
        trait_category3,
        perl = TRUE
      )
    ) |>
    mutate(
      trait_category3 = gsub(
        "(.*?):(.*?)",
        "\\1 - \\2",
        trait_category3,
        perl = TRUE
      )
    ) |>
    mutate(
      trait_category3 = str_remove(
        trait_category3,
        " ?two.component systems ?,?"
      )
    ) |>
    mutate(
      trait_category3 = str_replace(
        trait_category3,
        "desiccation/osmotic/salt stress",
        "osmotic stress"
      )
    ) |>
    mutate(
      trait_category3 = str_remove(
        trait_category3,
        "simple compound degradation - "
      )
    ) |>
    mutate(trait_category3 = str_remove(trait_category3, "C1 compounds - ")) |>
    mutate(trait_category3 = str_remove(trait_category3, "pigments - ")) |>
    mutate(
      trait_category3 = str_remove(
        trait_category3,
        "scavenging of reactive oxygen species - "
      )
    ) |>
    mutate(trait_category3 = str_remove(trait_category3, ".*transport - ")) |>
    mutate(
      trait_category3 = str_replace(
        trait_category3,
        "temperature - RNA",
        "Low temperature - RNA"
      )
    ) |>
    mutate(
      trait_category3 = str_remove(trait_category3, "chemo.*trophy - +")
    ) |>
    mutate(
      trait_category3 = if_else(
        trait_category3 %in%
          c(
            "ED pathway",
            "EMP pathway",
            "ETC complex IV",
            "Low temperature - RNA degradation",
            "OH-spheroidenone",
            "pH stress - redox sensing",
            "pH stress - urease system",
            "ROS scavenging enzymes",
            "S compound transport",
            "TCA cycle"
          ),
        trait_category3,
        str_to_sentence(trait_category3)
      )
    ) |>
    mutate(
      cat_facet = case_when(
        trait_category1 == "Stress Tolerance" ~ trait_category1,
        .default = trait_category2
      )
    ) |>
    mutate(cat_facet = str_to_sentence(cat_facet))

  data_clust <- plot_data |>
    select(env_label_good, trait_category3, per_genomes) |>
    pivot_wider(values_from = per_genomes, names_from = trait_category3) |>
    column_to_rownames("env_label_good")

  data.dist <- vegdist(t(data_clust), method = "euclidian", na.rm = TRUE)
  col.clus <- hclust(data.dist, "aver")
  ord <- col.clus$order

  order_cat <- colnames(data_clust)[ord]

  plot_data$trait_category3 <- factor(
    plot_data$trait_category3,
    levels = order_cat
  )

  heatmap_enriched_traits_full <- ggplot(
    plot_data,
    aes(x = env_label_good, y = trait_category3, fill = per_genomes)
  ) +
    geom_tile() +
    facet_grid(cat_facet ~ ., scales = "free", space = "free", switch = "y") +
    scale_fill_gradient2(
      midpoint = 50,
      low = "yellow",
      high = "red",
      mid = "orange",
      name = "Percentage of\ngenomes (%)"
    ) +
    theme_all +
    theme(
      axis.title = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(
        angle = 45,
        vjust = 1,
        hjust = 1
      )
    )

  heatmap_enriched_traits_full
}

# Joins a color palette to a node dataframe by environment label.
# Used to prepare collapse colors for the phylogenetic tree plot.
collapse_node_color <- function(df, color) {
  colors_df <- as.data.frame(color) |>
    rownames_to_column(var = "env_label_good")

  df_out <- df |>
    left_join(colors_df, join_by(env_label_good))

  df_out
}

# Draws the annotated phylogenetic tree with ancestral state colors.
# Collapses selected clades and overlays posterior probability sizes.
phylogenetic_tree <- function(
  rooted_tree,
  nodes_interest,
  max_prob,
  colors_samples,
  collapse_color
) {
  nodes_plot <- nodes_interest |>
    rename(nodes = node) |>
    mutate(reason = str_to_title(reason))

  colors_mod <- colors_samples

  names(colors_mod)[
    names(colors_mod) == "Other"
  ] <- "Uncertain"

  colors_mod[colors_mod == "black"] <- "gray44"

  max_prob_cat <- max_prob |>
    ungroup() |>
    mutate(
      cat = case_when(
        value >= 0.9 ~ "> 90%",
        value < 0.7 ~ "< 70%",
        .default = "70-90%"
      )
    ) |>
    select(node, cat)

  rooted_tree@data <- rooted_tree@data |>
    mutate(
      shape = if_else(node == 283, "Root", "Normal")
    ) |>
    mutate(
      env_label_good = if_else(
        env_label_good == "Lakewater",
        "Lake water",
        env_label_good
      )
    ) |>
    mutate(
      env_label_good = if_else(
        env_label_good == "Other",
        "Uncertain",
        env_label_good
      )
    ) |>
    left_join(
      max_prob_cat,
      join_by(node)
    )

  size_value = c(3, 4, 5)
  names(size_value) <- c("< 70%", "70-90%", "> 90%")

  p1 <- ggtree(rooted_tree, size = 0.3, alpha = 0.7) +
    geom_rootedge(rootedge = 10) +
    geom_nodepoint(aes(color = env_label_good, size = cat)) +
    # scale_shape_manual(values = c(16, 13), guide = "none") +
    scale_size_manual(
      values = size_value,
      name = "Posterior probabilities",
      breaks = c("< 70%", "70-90%", "> 90%")
    ) +
    scale_color_manual(
      values = colors_mod,
      na.value = "gray44",
      guide = "none"
    ) +
    geom_fruit(
      geom = geom_tile,
      aes(fill = env_label_good),
      width = 20,
      offset = 0.05
    ) +
    scale_fill_manual(values = colors_mod, name = "Environment") +
    new_scale_fill() +
    new_scale_color()

  p2 <- p1 +
    geom_point2(
      mapping = aes(subset = (node %in% nodes_plot$nodes), size = cat),
      shape = 21,
      # size = 5,
      stroke = 1.5
    ) +
    scale_fill_manual(values = colors_mod, name = "Environment")

  p3 <- reduce2(
    collapse_color$node,
    collapse_color$color,
    .init = p2,
    .f = function(plot, node, color) {
      collapse(plot, node, "mixed", fill = color, alpha = 0.8)
    }
  )
  p3
}

# Creates a grouped bar chart of DTL event counts at key nodes.
# Compares selected nodes against the average background.
plot_info_dtl <- function(DTL_nodes_int, colors_DTL) {
  p1 <- DTL_nodes_int |>
    pivot_longer(cols = !c(reason, node)) |>
    arrange(node) |>
    filter(!(reason %in% c("max", "min"))) |>
    mutate(
      reason = case_when(
        reason == "mean" ~ "Average",
        reason == "GFS" ~ "GFS",
        reason == "Glacier" ~ "Soil -> Glacier",
        reason == "GroundWater" ~ "Groundwater",
        reason == "LakeWater" ~ "Lake water",
        reason == "Root" ~ "LCA",
        reason == "Soil" ~ "Soil",
        reason == "Wetland" ~ "Lake water -> Wetland"
      )
    ) |>
    mutate(
      reason = factor(
        reason,
        levels = c(
          "Average",
          "LCA",
          "GFS",
          "Groundwater",
          "Lake water",
          "Lake water -> Wetland",
          "Soil",
          "Soil -> Glacier"
        )
      )
    ) |>
    mutate(name = str_to_title(name)) |>
    mutate(
      name = factor(
        name,
        levels = c("Duplications", "Transfers", "Origination", "Losses")
      )
    ) |>
    ggplot(aes(x = reason, y = value, fill = name)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = colors_DTL) +
    theme_classic() +
    labs(x = "Nodes", y = "Number of Events", fill = "Type of Events") +
    theme_all +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_text(
        angle = 45,
        vjust = 1,
        hjust = 1
      )
    )

  p1
}

# Plots gene prevalence distributions per pangenome partition category.
# Shows the percentage of genomes carrying each gene family.
plot_prevalence_ppan <- function(
  gene_pres_abs,
  partitions
) {
  dat <- gene_pres_abs |>
    pivot_longer(cols = !"Gene") |>
    group_by(Gene) |>
    summarise(n = sum(value), perc = n * 100 / 282) |>
    ungroup() |>
    left_join(partitions, join_by(Gene == Family)) |>
    mutate(Category = str_to_title(Category))

  ggplot(dat, aes(x = Category, y = perc, fill = Category)) +
    geom_boxplot() +
    labs(
      y = "Prevalence in\nPolaromonas genones (%)",
      x = "Type",
      fill = "Type"
    ) +
    theme_classic() +
    theme_all
}

# Draws a faceted heatmap of DTL events linked to microtrait functions.
# Rows are trait types, columns are event types, faceted by category and node.
heatmap_micro_dtl <- function(ev_micro_ancester_sel) {
  dat <- ev_micro_ancester_sel |>
    filter(presenceRemovingLosses > 0) |>
    select(
      microtrait_rule_name,
      "microtrait_trait_name3",
      duplications:losses,
      presenceRemovingLosses:reason
    ) |>
    distinct() |>
    select(-c(microtrait_rule_name)) |>
    pivot_longer(cols = c(duplications:presenceRemovingLosses)) |>
    group_by(across(-value)) |>
    summarise(sum = sum(value)) |>
    ungroup() |>
    filter(
      microtrait_trait_name3 !=
        "Resource Acquisition:Substrate degradation:simple compound degradation:protein degradation"
    ) |>
    separate_wider_delim(
      microtrait_trait_name3,
      delim = ":",
      names = c("c1", "c2", "c3", "c4", "c5", "c6", "c7"),
      too_few = "align_start",
      cols_remove = FALSE
    ) |>
    mutate(
      type = case_when(
        # !is.na(microtrait_rule_substrate) ~ paste(c3, microtrait_rule_substrate, sep = " - "),
        is.na(c4) ~ c3,
        is.na(c5) ~ c4,
        is.na(c6) ~ paste(c3, c5, sep = " - "),
        c3 == "temperature" ~ paste(c4, " temperature - ", c6, sep = ""),
        c3 == "pigments" ~ paste(c3, c6, sep = " - "),
        .default = c4
      )
    ) |>
    mutate(
      type = str_replace(
        type,
        "desiccation/osmotic/salt stress",
        "Osmotic stress"
      )
    ) |>
    mutate(type = str_remove(type, "simple compound degradation - ")) |>
    mutate(type = str_remove(type, "C1 compounds - ")) |>
    mutate(type = str_remove(type, "pigments - ")) |>
    mutate(type = str_remove(type, ".*transport - ")) |>
    mutate(type = str_remove(type, "chemo.*trophy - +")) |>
    mutate(
      category = case_when(
        c1 == "Resource Acquisition" ~ c2,
        c1 == "Resource Use" ~ c2,
        c1 == "Stress Tolerance" ~ "Stress Tolerance"
      )
    ) |>
    complete(type, reason) |>
    filter(!is.na(c1)) |>
    group_by(category, type, reason) |>
    mutate(nb = sum(sum)) |>
    ungroup() |>
    filter(nb > 0) |>
    mutate(sum = if_else(sum > 0, sum, NA)) |>
    # filter(reason != "Root") |>
    mutate(
      name = if_else(name == "presenceRemovingLosses", "Genes number", name)
    ) |>
    mutate(name = str_to_sentence(name)) |>
    mutate(
      name = factor(
        name,
        levels = c(
          "Genes number",
          "Duplications",
          "Transfers",
          "Origination",
          "Losses"
        )
      )
    ) |>
    mutate(
      reason = case_when(
        reason == "GFS" ~ "GFS",
        reason == "Glacier" ~ "Soil -> Glacier",
        reason == "GroundWater" ~ "Groundwater",
        reason == "LakeWater" ~ "Lake water",
        reason == "Root" ~ "LCA",
        reason == "Soil" ~ "Soil",
        reason == "Wetland" ~ "Lake water -> Wetland"
      )
    ) |>
    mutate(
      reason = factor(
        reason,
        levels = c(
          # "LUCA",
          "LCA",
          "GFS",
          "Groundwater",
          "Lake water",
          "Lake water -> Wetland",
          "Soil",
          "Soil -> Glacier"
        )
      )
    ) |>
    mutate(
      type = case_when(
        str_detect(type, "high temperature") ~
          str_replace(type, "high", "High"),
        str_detect(type, "low temperature") ~ str_replace(type, "low", "Low"),
        str_detect(type, "Osmotic") ~ type,
        str_detect(type, "pH") ~ type,
        str_detect(type, "photosytem") ~ str_replace(type, "photo", "Photo"),
        str_detect(type, "^oxidative") ~ str_replace(type, "oxi", "Oxi"),
        str_detect(type, "vitam") ~ str_replace(type, "vita", "Vita"),
        type %in%
          c(
            "ATP-dependent proteases",
            "ED pathway",
            "EPS biosynthesis/export",
            "PHB cycle",
            "ROS scavenging enzymes",
            "S compound transport"
          ) ~
          type,
        .default = str_to_sentence(type)
      )
    ) |>
    mutate(category = str_to_sentence(category)) |>
    mutate(sum = if_else(is.na(sum), 0, sum)) |>
    group_by(reason, category, type, name) |>
    summarise(temp = sum(sum)) |>
    ungroup() |>
    complete(reason, nesting(category, type), name, fill = list(temp = 0)) |>
    mutate(sum = if_else(temp == 0, NA, temp)) |>
    select(-temp) |>
    ungroup()

  p1 <- ggplot(
    dat,
    aes(
      x = reason,
      y = fct_rev(type),
      fill = sum
    )
  ) +
    geom_tile(color = "black") +
    facet_grid(
      category ~ name,
      scales = "free",
      space = "free",
      switch = "y"
    ) +
    scale_fill_gradient(
      low = "#fed976",
      high = "#7f0000",
      na.value = "white",
      name = "Number of Genes"
    ) +
    theme_all +
    theme(
      strip.text.x.top = element_text(angle = 0),
      axis.text.x = element_text(
        angle = 45,
        vjust = 1,
        hjust = 1
      ),
      axis.title = element_blank()
    )

  p1
}

# Plots the percentage of genomes carrying virus or plasmid signals per environment.
# Uses a dodged bar chart ordered by environment type.
plot_genomad_perc <- function(genomad_numb) {
  genomad_numb <- genomad_numb |>
    mutate(
      env_n = factor(
        env_n,
        levels = c(
          "Glacier (73)",
          "GFS (54)",
          "River (17)",
          "Wetland (26)",
          "Lake water (55)",
          "Groundwater (23)",
          "Soil (31)",
          "Other (3)"
        )
      )
    )

  ggplot(genomad_numb, aes(x = env_n, y = n, fill = type)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    theme_classic() +
    labs(y = "Percentage of\nGenomes (%)", fill = "Extrachromosomal DNA") +
    theme_all +
    theme(axis.title.x = element_blank())
}

# Plots DTL event frequencies linked to genomad elements at transition nodes.
# Facets by element type (virus/plasmid) with stacked bars per event type.
plot_genomad_dtl <- function(dat, colors) {
  dat_plot <- dat |>
    mutate(perc = perc * 100) |>
    mutate(type = str_to_title(type)) |>
    mutate(
      name = case_when(
        name == "dup_sum" ~ "Duplications",
        name == "orig_sum" ~ "Origination",
        name == "loss_sum" ~ "Losses",
        name == "trans_sum" ~ "Transfers"
      )
    ) |>
    mutate(
      name = factor(
        name,
        levels = c("Duplications", "Transfers", "Origination", "Losses")
      )
    ) |>
    mutate(
      reason = case_when(
        reason == "GFS" ~ "GFS",
        reason == "Glacier" ~ "Soil -> Glacier",
        reason == "GroundWater" ~ "Groundwater",
        reason == "LakeWater" ~ "Lake water",
        reason == "Root" ~ "LCA",
        reason == "Soil" ~ "Soil",
        reason == "Wetland" ~ "Lake water -> Wetland"
        # reason == "GFS" ~ reason,
        # reason == "LakeWater" ~ "Lake water",
        # .default = str_to_title(reason)
      )
    ) |>
    mutate(reason_numb = paste(reason, "\n(", allnumb, ")", sep = "")) |>
    mutate(
      reason_numb = factor(
        reason_numb,
        levels = c(
          "LCA\n(23)",
          "LCA\n(820)",
          "GFS\n(2013)",
          "GFS\n(59)",
          "Groundwater\n(1891)",
          "Groundwater\n(50)",
          "Lake water\n(1842)",
          "Lake water\n(52)",
          "Lake water -> Wetland\n(1845)",
          "Lake water -> Wetland\n(54)",
          "Soil\n(2220)",
          "Soil\n(63)",
          "Soil -> Glacier\n(2334)",
          "Soil -> Glacier\n(70)"
        )
      )
    )

  p1 <- ggplot(dat_plot, aes(x = reason_numb, y = perc, fill = name)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = colors) +
    facet_wrap(~type, scales = "free_x") +
    labs(y = "Frequence (%)", x = "", fill = "Type of Events") +
    theme_classic() +
    theme_all +
    theme(
      strip.text.x.top = element_text(angle = 0),
      axis.text.x = element_text(
        angle = 45,
        vjust = 1,
        hjust = 1
      ),
    )

  p1
}

# Plots a rarefaction curve of pangenome gene categories vs. genome count.
# Uses log10 y-axis with median and SD error bars.
rarefaction_curve <- function(dat) {
  rare_plot <- dat |>
    select(genomes_count:cloud) |>
    na.omit() |>
    pivot_longer(cols = !genomes_count) |>
    mutate(name = str_to_title(name)) |>
    filter(value > 0) |>
    ggplot(aes(x = genomes_count, y = value, color = name)) +
    # geom_point()#+
    stat_summary(fun = median, geom = "point", size = 1, shape = 18) +
    scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
    stat_summary(
      fun.data = function(x) {
        data.frame(ymin = median(x) - sd(x), ymax = median(x) + sd(x))
      },
      geom = "errorbar",
      width = 0.2,
      size = 0.5
    ) +
    labs(
      x = "Number of Genomes",
      y = "Number of Genes",
      color = "Genes Category"
    ) +
    theme_classic() +
    theme_all

  rare_plot
}

# Draws a circular phylogenetic tree colored by neighboring taxon groups.
# Adds outer rings for taxon identity and monophyly status.
plot_neigh_monophy <- function(tree, mono) {
  getPaletteBact = colorRampPalette(brewer.pal(12, "Paired"))
  treeColor <- getPaletteBact(length(unique(mono$Taxon)) + 1)
  treeColor[1] <- "black"

  groups <- split(mono$Tip, mono$Taxon)
  tree_groups <- groupOTU(tree, groups)

  p1 <-
    ggtree(tree_groups, layout = 'circular', aes(color = group)) + #
    geom_tree() +
    theme_tree() +
    scale_color_manual(
      values = treeColor,
      na.value = "transparent",
      guide = "none"
    ) +
    geom_fruit(
      data = mono,
      pwidth = 0.03,
      offset = 0.05,
      geom = geom_bar,
      mapping = aes(y = Tip, fill = Taxon, x = 1),
      stat = "identity",
    ) +
    scale_fill_manual(values = treeColor[-1]) +
    new_scale_colour() +
    new_scale_fill() +
    geom_fruit(
      data = mono,
      pwidth = 0.03,
      offset = 0.02,
      geom = geom_bar,
      mapping = aes(y = Tip, fill = Status, x = 1),
      stat = "identity",
    )

  p1
}

# Formats the variance partitioning results into a gt table.
# Displays adjusted R² and p-values for environment and phylogeny fractions.
format_gene_phylo <- function(df) {
  df_gt <- df |>
    group_by(type) |>
    gt(rowname_col = "variable") |>
    fmt_percent(Adj.R.squared, decimals = 1) |>
    sub_missing(missing_text = "") |>
    sub_zero(columns = Df, zero_text = "") |>
    cols_label(
      "Df" = "Degree of\nFreedom",
      "Adj.R.squared" = "Adjusted {{R^2}}",
      "p_value" = "p-value"
    ) |>
    cols_align(align = "center", columns = !"variable")

  df_gt
}

# Formats PERMANOVA results into a gt table grouped by partition type.
# Shows Df, R², F-statistic, and p-value per model term.
format_gene_complete <- function(df) {
  df_gt <- df |>
    select(-SumOfSqs) |>
    filter(variable != "Total") |>
    group_by(type) |>
    gt(rowname_col = "variable") |>
    fmt_percent(R2, decimals = 1) |>
    fmt_number("F", decimals = 1) |>
    sub_missing(missing_text = "") |>
    sub_zero(columns = Df, zero_text = "") |>
    cols_label(
      "Df" = "Degree of\nFreedom",
      "R2" = "{{R^2}}",
      "Pr(>F)" = "p-value"
    ) |>
    cols_align(align = "center", columns = !"variable")

  df_gt
}

# Compares node-level confidence across ACE, pastML, and stochastic mapping.
# Plots histograms of max marginal probabilities with count annotations.
plot_node_confidence <- function(ace, pastml, stoch) {
  df_ace <- ace |>
    select(label, Max_prob, Top_state) |>
    rename(value = Max_prob) |>
    mutate(type = "Ancestral Character Estimation")

  df_past <- pastml |>
    select(label, name, value) |>
    rename(Top_state = name) |>
    mutate(type = "pastML")

  df_stoch <- stoch |>
    select(ML_label, SIMMAP_state, SIMMAP_pp) |>
    rename(c(Top_state = SIMMAP_state, value = SIMMAP_pp, label = ML_label)) |>
    mutate(type = "Stochastic mapping")

  df <- rbind(df_ace, df_past, df_stoch) |>
    mutate(value = value * 100)

  df_summary <- df |>
    group_by(type) |>
    summarize(
      less_7 = sum(value < 70),
      above_90 = sum(value >= 90),
      other = 281 - less_7 - above_90
    ) |>
    ungroup() |>
    pivot_longer(cols = !type) |>
    mutate(
      x_value = case_when(
        name == "less_7" ~ 55,
        name == "above_90" ~ 95,
        .default = 80
      )
    )

  p <- ggplot() +
    geom_histogram(data = df, aes(value, fill = type)) +
    geom_text(
      data = df_summary,
      aes(label = value, x = x_value),
      y = 150,
      size = 4.5
    ) +
    geom_vline(xintercept = 70, linetype = 3, col = "red", linewidth = 0.7) +
    geom_vline(
      xintercept = 90,
      linetype = 3,
      col = "darkgreen",
      linewidth = 0.7
    ) +
    facet_wrap(~type) +
    guides(fill = "none") +
    labs(
      x = "Maximum marginal probability (%)",
      y = "Number of nodes",
      title = "Ancestral node confidence (ML marginal posteriors)",
    ) +
    theme_classic() +
    theme_all +
    theme(
      strip.text.x.top = element_text(angle = 0)
    )
  p
}

# Formats enrichment test results for MGE transfers into a gt table.
# Reports permutation p-value, effect size, and bootstrap confidence interval.
format_table_nodes_transfers <- function(dat) {
  dat_fmt <- dat |>
    arrange(label) |>
    mutate(
      label = case_when(
        str_detect(label, "plasmid") ~ "MGE - Plasmid",
        str_detect(label, "all") ~ "MGE - All",
        str_detect(label, "virus") ~ "MGE - Virus",
        .default = label
      )
    ) |>
    gt() |>
    fmt_scientific(perm) |>
    fmt_number(columns = c(eff_size, ci_low, ci_high), decimals = 2) |>
    tab_spanner(
      label = "CI",
      columns = c(ci_low, ci_high)
    ) |>
    cols_align(align = "center", columns = !label) |>
    cols_label(
      perm = "Permutation p",
      eff_size = "Effect size (r)",
      ci_low = "Low",
      ci_high = "High",
      label = ""
    )

  dat_fmt
}

# Formats significant PGLS (caper) results into a styled gt table.
# Cleans trait names and displays lambda, LRT, and p-value per trait.
table_caper <- function(dat) {
  dat_fmt <- dat |>
    na.omit() |>
    filter(p_value <= 0.05) |>
    separate_wider_delim(
      pathway,
      delim = "__",
      names = c("c1", "c2", "c3", "c4", "c5", "c6", "c7"),
      too_few = "align_start",
      cols_remove = TRUE
    ) |>
    mutate(
      type = case_when(
        is.na(c4) ~ c3,
        is.na(c5) ~ c4,
        is.na(c6) ~ paste(c3, c5, sep = " - "),
        c3 == "temperature" ~ paste(c4, " temperature - ", c6, sep = ""),
        c3 == "pigments" ~ paste(c3, c6, sep = " - "),
        .default = c4
      )
    ) |>
    mutate(type = str_remove(type, "simple_compound_degradation - ")) |>
    mutate(type = str_remove(type, "C1_compounds - ")) |>
    mutate(type = str_remove(type, "pigments - ")) |>
    mutate(type = str_remove(type, ".*transport - ")) |>
    mutate(type = str_remove(type, "chemo.*trophy - +")) |>
    mutate(
      category = case_when(
        c1 == "Resource_Acquisition" ~ c2,
        c1 == "Resource_Use" ~ c2,
        c1 == "Stress_Tolerance" ~ "Stress_Tolerance"
      )
    ) |>
    select(c(category, type, lambda, LRT_stat, p_value, signal_class)) |>
    arrange(category, type) |>
    gt() |>
    text_transform(
      fn = function(x) {
        gsub("_", " ", x)
      },
      locations = cells_body()
    ) |>
    text_transform(
      fn = function(x) {
        stringr::str_to_sentence(x)
      },
      locations = cells_body()
    ) |>
    fmt_scientific(p_value) |>
    fmt_number(columns = c(lambda, LRT_stat), decimals = 2) |>
    cols_label(
      category = "Category",
      type = "Type",
      lambda = paste("Pagel's", "\U03BB"),
      LRT_stat = "Likelihood ratio test",
      p_value = "p value",
      signal_class = "Meaning"
    ) |>
    cols_align(align = "center", columns = !c("category", "type"))

  dat_fmt
}

# Formats the sensitivity analysis results of ACE into a gt table.
# Reports mean and SD agreement across perturbation rates.
table_sensitivity <- function(dat) {
  dat_fmt <- dat |>
    select(-Min_agreement) |>
    gt() |>
    fmt_percent(columns = !Rate, decimals = 1) |>
    fmt_percent(columns = Rate, decimals = 0) |>
    tab_spanner(
      label = "Agreement",
      columns = c(Mean_agreement, SD_agreement)
    ) |>
    cols_label(
      Rate = "Misclassification\nrate",
      Mean_agreement = "Average",
      SD_agreement = "Std"
    ) |>
    cols_align(align = "center")

  dat_fmt
}
