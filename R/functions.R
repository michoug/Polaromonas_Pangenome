# Cleans and processes environmental metadata

clean_env_metadata <- function(dat) {
  dat_clean <- dat |>
    mutate(.by = env_label_manual_2, n = n_distinct(user_genome)) |>
    mutate(
      env_label_good = case_when(
        is.na(env_label_manual_2) ~ "Other",
        env_label_manual_2 == "lakewater" ~ "Lake water",
        env_label_manual_2 == "mine drainage" ~ "Other",
        env_label_manual_2 == "Saline lake" ~ "Lake water",
        env_label_manual_2 == "glacier sediment" ~ "Glacier",
        env_label_manual_2 == "ice_core" ~ "Glacier",
        env_label_manual_2 == "GFS" ~ "GFS",
        .default = str_to_title(env_label_manual_2)
      )
    ) |>
    select(user_genome, country, env_label_good) |>
    mutate(
      env_label_good = factor(
        env_label_good,
        levels = c(
          "Glacier",
          "GFS",
          "River",
          "Wetland",
          "Lake water",
          "Groundwater",
          "Soil",
          "Other"
        )
      )
    ) |>
    rename("Genome" = "user_genome") |>
    group_by(country) |>
    mutate(n = n()) |>
    ungroup() |>
    mutate(country_good = if_else(n > 5, country, "Other"))
}

combine_metadata <- function(meta, stats, check) {
  dat_clean <- meta |>
    mutate(
      Environment = case_when(
        is.na(env_label_manual_2) ~ "Other",
        env_label_manual_2 == "lakewater" ~ "Lake water",
        env_label_manual_2 == "mine drainage" ~ "Other",
        env_label_manual_2 == "Saline lake" ~ "Lake water",
        env_label_manual_2 == "glacier sediment" ~ "Glacier",
        env_label_manual_2 == "ice_core" ~ "Glacier",
        env_label_manual_2 == "GFS" ~ "GFS",
        .default = str_to_title(env_label_manual_2)
      )
    ) |>
    left_join(check, join_by("user_genome" == "Name")) |>
    rename("Genome" = "user_genome") |>
    left_join(stats, join_by("Genome")) |>
    select(
      -c(env_label_manual_2, "Completeness_Model_Used", "Additional_Notes")
    ) |>
    rename(c(
      "Assembly accession" = "assembly_accession",
      "Bioproject" = "assembly_bio_project_accession",
      "Contigs" = "# Sequences",
      "Country" = "country"
    ))
}

# Always use the same colors for the environment!
getColors <- function() {
  colors <- c(brewer.pal(7, "Accent"), "black")
  names(colors) <- c(
    "Glacier",
    "GFS",
    "River",
    "Wetland",
    "Lake water",
    "Groundwater",
    "Soil",
    "Other"
  )
  # colors <-
  #   c(
  #     "Glacier" = "#1F78B4",
  #     "GFS" = "#CAB2D6",
  #     "River" = "#6A3D9A",
  #     "Wetland" = "#B2DF8A",
  #     "Lake water" = "#33A02C",
  #     "Groundwater" = "#FF7F00",
  #     "Soil" = "#B15928",
  #     "Other" = "darkgrey"
  #   )

  colors[4] <- "#ffee99"
  colors
}

getColorsCountry <- function(dat) {
  colors_country <- c(
    rgb(235, 172, 35, maxColorValue = 255),
    rgb(184, 0, 88, maxColorValue = 255),
    rgb(0, 140, 249, maxColorValue = 255),
    rgb(0, 110, 0, maxColorValue = 255),
    rgb(0, 187, 173, maxColorValue = 255),
    rgb(209, 99, 230, maxColorValue = 255),
    rgb(178, 69, 2, maxColorValue = 255),
    rgb(255, 146, 135, maxColorValue = 255),
    rgb(89, 84, 214, maxColorValue = 255),
    rgb(0, 198, 248, maxColorValue = 255),
    "black"
  )

  country_name <- unique(dat$country_good) |>
    sort()

  country_name <- c(setdiff(country_name, "Other"), "Other")

  names(colors_country) <- country_name

  colors_country
}

getcolorsDTL <- function() {
  colors_dtl <- c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3")

  names(colors_dtl) <- c(
    "Duplications",
    "Transfers",
    "Origination",
    "Losses"
  )
  colors_dtl
}


combine_atlas_df <- function(zipfile, folder) {
  if (!dir.exists(folder)) {
    dir.create(folder)
  }

  unzip(zipfile, exdir = folder, overwrite = TRUE)

  lst <- dir_ls(folder) |>
    set_names(basename) |>
    purrr::map(read_tsv, show_col_types = FALSE) |>
    list_rbind(names_to = "type")

  unlink(folder, recursive = TRUE, force = TRUE)
  lst
}

atlas_DF_map <- function(df, coord, samples) {
  df_clean <- df |>
    clean_names() |>
    mutate(samples = paste(run_id, number_sample_id, sep = ".")) |>
    left_join(samples, join_by("samples")) |>
    mutate(
      nbRead = abundance_ppm * RNA_reads / 1e6,
      percAbund = nbRead * 100 / RNA_reads
    ) |>
    filter(percAbund > 0.1) |>
    mutate(
      env_good = case_when(
        str_detect(environment, "animal") ~ "animal",
        str_detect(environment, "ice") ~ "ice",
        str_detect(environment, "river") ~ "river",
        str_detect(environment, "groundwater") ~ "groundwater",
        str_detect(environment, "lake") ~ "lake",
        str_detect(environment, "sediment") ~ "sediment",
        str_detect(environment, "soil") ~ "soil",
        str_detect(environment, "sea") ~ "marine",
        str_detect(environment, "marine") ~ "marine",
        str_detect(environment, "plant") ~ "plant",
        str_detect(environment, "waste water") ~ "waste water",
        str_detect(keywords, "glacial") ~ "glacier",
        is.na(environment) ~ "Unknown",
        .default = environment
      )
    ) |>
    mutate(env_good = str_to_sentence(env_good)) |>
    group_by(number_sample_id, samples, env_good) |>
    summarise(percAbund = sum(percAbund)) |>
    ungroup() |>
    left_join(coord, join_by(number_sample_id == sampleId)) |>
    filter(parsedLon != "None") |>
    select(
      number_sample_id,
      samples,
      percAbund,
      env_good,
      parsedLat,
      parsedLon
    ) |>
    filter(number_sample_id != "SRS7752820") |>
    mutate(parsedLon = as.numeric(parsedLon)) |>
    mutate(parsedLat = as.numeric(parsedLat)) |>
    filter(parsedLat < 90) |>
    mutate(env_good = as.factor(env_good)) |>
    mutate(env_good = factor(env_good, levels = sort(unique(env_good)))) |>
    mutate(env_good = fct_relevel(env_good, "Unknown", after = Inf))

  df_clean
}

combine_stats_checkm_temp <- function(stats, check, temp) {
  dat <- stats |>
    left_join(check, join_by(Genome == Name)) |>
    left_join(temp, join_by(Genome)) |>
    mutate(relative_length = Size * (100 / Completeness))
}

combine_with_groups <- function(dat, cat) {
  dat_g <- dat |>
    left_join(cat, join_by("Genome"))
}

get_contig_prot <- function(zipfile, folder) {
  if (!dir.exists(folder)) {
    dir.create(folder)
  }

  unzip(zipfile, exdir = folder, overwrite = TRUE)

  lst <- dir_ls(folder) |>
    set_names(basename) |>
    purrr::map(
      read_tsv,
      show_col_types = FALSE,
      progress = FALSE,
      col_types = "cc",
      comment = "#",
      col_names = c("Contig", "Protein")
    ) |>
    list_rbind(names_to = "Genome") |>
    mutate(
      Genome = Genome |>
        str_replace_all(".tsv", "")
    )

  unlink(folder, recursive = TRUE, force = TRUE)
  lst
}

get_genomad_dtl <- function(prot, families, ev, virus, plasmid) {
  vir <- virus |>
    mutate(type = "virus") |>
    select(Genome, seq_name, type)

  plas <- plasmid |>
    mutate(type = "plasmid") |>
    select(Genome, seq_name, type)

  dat <- vir |>
    rbind(plas) |>
    select(-Genome) |>
    left_join(
      prot,
      join_by("seq_name" == "Contig"),
      relationship = "many-to-many"
    ) |>
    left_join(families, join_by("Protein" == "Gene")) |>
    select(-c("localGene", "Frag", "Protein")) |>
    distinct() |>
    na.omit() |>
    left_join(ev, join_by("Family"), relationship = "many-to-many")
}

prepare_genes_for_nmds <- function(genes, partition, category) {
  partition_select <- partition |>
    filter(Category == category)

  dat <- genes |>
    filter(Gene %in% partition_select$Family) |>
    column_to_rownames("Gene") |>
    select(where(function(x) sum(x) > 0)) |>
    as.matrix() |>
    t()

  dat
}

create_nmds <- function(dat) {
  nmds <- metaMDS(
    dat,
    k = 2,
    trace = 1,
    autotransform = FALSE,
    distance = "jaccard",
    trymax = 100
  )

  nmds
}

anosim_nmds <- function(dat, map) {
  rownames(dat) <- gsub("-", "_", rownames(dat))

  map_clean <- map |>
    arrange(match(Genome, rownames(dat)))

  dat_dist <- (dat) |>
    vegdist(method = "jaccard")

  dat.ano <- with(
    map_clean,
    anosim(dat_dist, env_label_good, distance = "jaccard")
  )

  dat.ano
}


clean_env_stats <- function(env, stats) {
  dat <- env |>
    left_join(stats, join_by("Genome")) |>
    mutate(Genome = gsub("-", "_", Genome)) |>
    select(Genome, env_label_good, relative_length, GC, CDS)
}

prepare_microtrait_data <- function(micro, stats) {
  microtrait_clean <- micro$trait_matrixatgranularity3 |>
    mutate(user_genome = gsub("_new", "", id), .keep = "unused") |>
    mutate(Genome = gsub("-", "_", user_genome), .keep = "unused") |>
    relocate(Genome) |>
    left_join(stats, join_by(Genome)) |>
    pivot_longer(
      cols = "Resource Acquisition:Substrate uptake:aromatic acid transport":"Stress Tolerance:Specific:envelope stress:phage resistance:phage-shock system"
    ) |>
    mutate(value = if_else(value >= 1, 1, 0)) |>
    mutate(name = gsub(" ", "_", name)) |>
    mutate(name = gsub(":", "__", name)) |>
    mutate(name = gsub("\\(|\\)|,", "", name))

  microtrait_clean
}

stats_microtrait_enrich <- function(micro) {
  env_gradient_df <- data.frame(
    env_label_good = c("Glacier", "GFS", "River", "Wetland", "Lake water"),
    env_gradient_num = c(1, 2, 3, 4, 5)
  )

  microtrait_stats <- micro |>
    left_join(env_gradient_df, join_by(env_label_good)) |>
    drop_na(env_gradient_num) |>
    group_by(name) |>
    filter(sum(value) > 0) |>
    summarise(
      p.val = summary(glm(
        value ~ CDS + env_gradient_num,
        family = binomial(link = "logit")
      ))$coefficients["env_gradient_num", 4] # test if environments has significant effect
      #,CDS_p.val = summary(glm(value~ CDS + env_gradient_num,family = "binomial"))$coefficients["CDS",4] # test if the number of CDSs has significant effect
    ) |>
    ungroup() |>
    mutate(p.adj = round(p.adjust(p.val, "BH"), 3)) |>
    filter(p.adj <= 0.05)

  microtrait_stats
}

combine_events_df <- function(zipfile, directory = tempfile()) {
  if (!dir.exists(directory)) {
    dir.create(directory)
  }

  unzip(zipfile, exdir = directory, overwrite = TRUE)

  indivEventlist <- fs::dir_ls(directory, glob = "*count.txt.gz")

  indivEvent <- vroom(
    indivEventlist,
    id = "family",
    show_col_types = FALSE,
    delim = "," # or "\t" if tab-separated
  )

  unlink(directory, recursive = TRUE, force = TRUE)

  return(indivEvent)
}

clean_egg <- function(dat) {
  eggnog_ko <- dat |>
    filter(KEGG_ko != "-") |>
    separate_rows(KEGG_ko, sep = ",") |>
    mutate(KEGG_ko = gsub("ko:", "", KEGG_ko)) |>
    select(query, KEGG_ko)
  eggnog_ko
}

get_families_Kos <- function(dat, egg) {
  dat_clean <- dat |>
    select(Family) |>
    distinct() |>
    inner_join(egg, join_by(Family == query))

  dat_clean
}


clean_events <- function(ev) {
  events_clean <- ev |>
    filter(str_detect(species_label, "Node")) |>
    filter(speciations > 0.3) |>
    mutate(Family = gsub("_perspecies_eventcount.txt.gz", "", family)) |>
    mutate(Family = gsub(".nex", "", Family)) |>
    mutate(Family = gsub("temp/", "", Family)) |>
    select(
      Family,
      species_label,
      duplications,
      transfers,
      origination,
      losses
    ) |>
    mutate(across(duplications:losses, ~ if_else(.x > 0.3, 1, 0)))

  events_clean
}


get_nodes_tree <- function(tree) {
  p1 <- ggtree(tree)

  node_data <- p1$data |>
    filter(!isTip) |>
    select(node, label) |>
    filter(str_detect(label, "Node"))
  node_data
}

getDTLNodes <- function(dat, nodes, nodes_imp) {
  dat_info <- dat |>
    left_join(nodes, join_by(species_label == label)) |>
    right_join(nodes_imp, join_by(node)) |>
    group_by(node, reason) |>
    summarise(across(duplications:losses, ~ sum(.x)))

  dat_full <- dat |>
    left_join(nodes, join_by(species_label == label)) |>
    group_by(node) |>
    summarise(across(duplications:losses, ~ sum(.x))) |>
    ungroup() |>
    summarise(across(
      duplications:losses,
      list(mean = mean, min = min, max = max)
    )) |>
    pivot_longer(
      cols = everything(),
      names_to = c(".value", "reason"),
      names_sep = "_"
    ) |>
    mutate(node = 1) |>
    rbind(dat_info)

  dat_full
}

events_kos <- function(ev, fam_kos, node) {
  ev_kos <- ev |>
    left_join(fam_kos, join_by(Family), relationship = "many-to-many") |>
    group_by(species_label) |>
    mutate(nbFamilies = n_distinct(Family), nbKos = n_distinct(KEGG_ko)) |>
    select(!Family) |>
    ungroup() |>
    group_by(species_label, KEGG_ko, nbFamilies, nbKos) |>
    summarise(
      across(duplications:losses, ~ sum(.x)),
      presenceWithLosses = n(),
      .groups = "drop_last"
    ) |>
    na.omit() |>
    left_join(node, join_by(species_label == label))

  ev_kos
}

process_rds <- function(path) {
  df <- readRDS(path)

  df_clean <- df$rules_asserted |>
    clean_names() |>
    mutate(name = gsub(".*/(.*).microtrait.rds", "\\1", path)) |>
    filter(microtrait_rule_asserted) |>
    select(name, microtrait_rule_name, microtrait_rule_boolean) |>
    mutate(
      parsed_components = map(microtrait_rule_boolean, parse_boolean_expression)
    ) |>
    unnest(parsed_components) |>
    select(-microtrait_rule_boolean) |>
    distinct() |>
    as.data.frame()

  df_genes <- df$genes_detected_table |>
    select(gene_name, hmm_name)

  df_dup <- df_clean |>
    filter(microtrait_rule_name != parsed_components)

  list_dup <- intersect(df_dup$microtrait_rule_name, df_dup$parsed_components)

  df_clean_sel <- df_clean |>
    select(-name)

  df_non_dup <- df_clean |>
    full_join(
      df_clean_sel,
      join_by(parsed_components == microtrait_rule_name)
    ) |>
    mutate(
      parsed_components = if_else(
        is.na(parsed_components.y),
        parsed_components,
        parsed_components.y
      )
    ) |>
    select(-parsed_components.y) |>
    right_join(df_genes, join_by(parsed_components == hmm_name))

  df_non_dup
}

parse_boolean_expression <- function(expression) {
  cleaned_expression <- gsub("[()']", "", expression)
  components <- unlist(strsplit(cleaned_expression, "\\s*\\|\\s*|\\s*&\\s*"))
  components <- trimws(components)
  return(components)
}

combine_rds_df <- function(zipfile, folder) {
  if (!dir.exists(folder)) {
    dir.create(folder)
  }

  unzip(zipfile, exdir = folder, overwrite = TRUE)

  lst <- dir_ls(folder) |>
    purrr::map(process_rds) |>
    list_rbind()

  unlink(folder, recursive = TRUE, force = TRUE)
  lst
}

cleanmicro <- function(dat) {
  dat_clean <- dat |>
    mutate(hmm_name = `microtrait_rule-boolean`) |>
    separate_rows(hmm_name, sep = "&") |>
    separate_rows(hmm_name, sep = "\\|") |>
    mutate(hmm_name = str_remove_all(hmm_name, "['()' ]"))
}

clean_rules <- function(rule2, rule) {
  rule2clean <- rule2 |>
    cleanmicro()

  ruleclean <- rule |>
    cleanmicro() |>
    select(
      -c(
        "microtrait_rule-booleanunwrapped",
        "microtrait_rule-namequoted"
      )
    )

  rule_file <- ruleclean |>
    left_join(
      rule2clean,
      join_by(
        `microtrait_rule-name`,
        `microtrait_rule-boolean`,
        hmm_name
      )
    ) |>
    distinct() |>
    clean_names()
}

subset_nodes <- function(ev, node, node_int) {
  dat <- ev |>
    left_join(node, join_by(species_label == label)) |>
    right_join(node_int, join_by(node))
}

get_micro_events <- function(ev, rds, rules, node) {
  rds_clean <- rules |>
    mutate(microtrait_rule_name = as.character(microtrait_rule_name)) |>
    right_join(
      rds,
      join_by(microtrait_rule_name),
      relationship = "many-to-many"
    ) |>
    filter(!is.na(microtrait_trait_name3)) |>
    select(-name) |>
    distinct()

  dat_clean <- ev |>
    right_join(
      rds_clean,
      join_by(Family == gene_name),
      relationship = "many-to-many"
    ) |>
    group_by(species_label) |>
    mutate(nbFamiles = n_distinct(Family), nbhmm = n_distinct(hmm_name)) |>
    ungroup() |>
    select(-c(Family, parsed_components)) |>
    group_by(species_label, across(microtrait_rule_name:nbhmm)) |>
    summarise(
      across(duplications:losses, ~ sum(.x)),
      presenceWithLosses = n(),
      .groups = "drop_last"
    ) |>
    left_join(node, join_by(species_label == label))

  dat_clean
}

clean_kegg_mod <- function(dat, prok) {
  dat_clean <- dat |>
    select(KEGG_ko, Subcategory, Module, ModuleDescription) |>
    filter(Module %in% prok$Module)
  dat_clean
}

clean_kegg_path <- function(dat) {
  dat_clean <- dat |>
    filter(
      !(Category %in%
        c(
          "09150 Organismal Systems",
          "09160 Human Diseases",
          "09180 Brite Hierarchies",
          "09190 Not Included in Pathway or Brite"
        ))
    ) |>
    filter(
      !(Subcategory %in%
        c(
          "09143 Cell growth and death",
          "09123 Folding, sorting and degradation"
        ))
    ) |>
    mutate(Subcategory = gsub("^\\d+ ", "", Subcategory)) |>
    rename(KEGG_ko = KO)
  dat_clean
}

getFunctionAncesterKO <- function(ev, path, mod) {
  kegg_path_mod <- path |>
    full_join(mod, join_by(Subcategory, KEGG_ko), relationship = "many-to-many")

  ev_ancest_ko <- ev |>
    left_join(
      kegg_path_mod,
      join_by("KEGG_ko"),
      relationship = "many-to-many"
    ) |>
    group_by(Category, Subcategory, node, Pathway) |>
    mutate(UniqKos = n()) |>
    ungroup() |>
    na.omit() |>
    select(-c(Category, Enzyme, `EC Number`)) |>
    distinct() |>
    mutate(presenceRemovingLosses = presenceWithLosses - losses)

  ev_ancest_ko
}

getFunctionAncesterMicro <- function(dat) {
  dat_clean <- dat |>
    group_by(node, microtrait_trait_name3) |>
    mutate(nbGene = n()) |>
    ungroup() |>
    filter(!is.na(microtrait_trait_name3)) |>
    mutate(presenceRemovingLosses = presenceWithLosses - losses)

  dat_clean
}

selectFunctionAncester <- function(dat, node) {
  dat_sel <- dat |>
    right_join(node, join_by(node))
  dat_sel
}


get_genomad_dat <- function(zipfile, folder, type) {
  if (!dir.exists(folder)) {
    dir.create(folder)
  }

  unzip(zipfile, exdir = folder, overwrite = TRUE)

  lst <- dir_ls(folder) |>
    set_names(basename) |>
    purrr::map(
      read_tsv,
      show_col_types = FALSE,
      progress = FALSE,
      col_types = type
    ) |>
    list_rbind(names_to = "Genome") |>
    mutate(
      Genome = Genome |>
        str_replace_all("-", "_") |>
        str_remove("_(virus|plasmid)_summary\\.tsv$")
    )

  unlink(folder, recursive = TRUE, force = TRUE)
  lst
}

genomad_bar <- function(meta, virus, plasmid) {
  dat_virus <- virus |>
    select(Genome) |>
    mutate(type = "Virus") |>
    distinct()

  dat_plasmid <- plasmid |>
    select(Genome) |>
    mutate(type = "Plasmid") |>
    distinct()

  meta <- meta |>
    group_by(env_label_good) |>
    mutate(n_env = n()) |>
    ungroup() |>
    select(Genome, env_label_good, n_env)

  dat_all <- dat_virus |>
    rbind(dat_plasmid) |>
    left_join(meta, join_by("Genome")) |>
    mutate(env_n = glue("{env_label_good} ({n_env})")) |>
    group_by(type, env_label_good, env_n) |>
    reframe(n = n() * 100 / n_env) |>
    distinct()

  dat_all
}

genomad_dtl_filter <- function(genomad, tree, nodes) {
  genomad_dtl_filter <- genomad |>
    left_join(tree, join_by("species_label" == "label")) |>
    right_join(nodes, join_by("node")) |>
    mutate(Genome = gsub("-", "_", Genome)) |>
    group_by(type, reason) |>
    mutate(seq_genome = paste(Genome, seq_name, sep = "_")) |>
    mutate(allnumb = n_distinct(seq_genome)) |>
    ungroup() |>
    select(
      seq_genome,
      type,
      Family,
      duplications:reason,
      allnumb
    ) |>
    distinct() |>
    group_by(seq_genome, reason) |>
    mutate(
      n = n(),
      dup_sum = sum(duplications) / n,
      trans_sum = sum(transfers) / n,
      orig_sum = sum(origination) / n,
      loss_sum = sum(losses) / n
    ) |>
    ungroup() |>
    select(seq_genome, type, reason, allnumb, n:loss_sum) |>
    distinct() |>
    rowwise() |>
    mutate(max_row = max(c_across(dup_sum:loss_sum))) |>
    filter(max_row > 0.2) |>
    select(-c(max_row, n)) |>
    pivot_longer(cols = !c(seq_genome, type, reason, allnumb)) |>
    filter(value > 0) |>
    group_by(seq_genome, type, reason, allnumb) |>
    slice_max(order_by = value, n = 1, with_ties = FALSE) |>
    ungroup() |>
    group_by(type, reason, allnumb, name) |>
    reframe(perc = n() / allnumb) |>
    distinct()
}


monophyletic <- function(tree, tax) {
  tax_clean <- tax |>
    as.data.frame()

  solution1 <- AssessMonophyly(tree, tax_clean)

  tips_mono <- solution1$genus$TipStates
  tips_mono
}

gene_jaccard <- function(matrix, partition, cat) {
  if (cat == "All") {
    part_clean <- partition
  } else {
    part_clean <- partition |>
      filter(Category == cat)
  }

  gene_matrix <- matrix |>
    filter(Gene %in% part_clean$Family) |>
    pivot_longer(cols = !Gene) |>
    mutate(name = gsub("-", "_", name)) |>
    pivot_wider(names_from = "Gene") |>
    as.data.frame() |>
    column_to_rownames("name")

  gene_dist <- vegdist(gene_matrix, method = "jaccard", binary = TRUE)
}

# adonis_gene_phylogeny <- function(jaccard, tree, meta, cat) {
#   genome_ids <- rownames(jaccard)
#
#   phylo_dist <- cophenetic.phylo(tree@phylo)
#   phylo_dist <- as.dist(phylo_dist[genome_ids, genome_ids])
#
#   phylo_pcoa <- cmdscale(phylo_dist, k = 10, eig = TRUE)
#
#   phylo_axes <- as.data.frame(phylo_pcoa$points[, 1:5])
#   colnames(phylo_axes) <- paste0("PCoA", 1:5)
#
#   env_df <- meta |>
#   select(Genome, env_label_good) |>
#     as.data.frame() |>
#     arrange(match(Genome, genome_ids)) |>
#     column_to_rownames("Genome")
#
#   covariates <- cbind(phylo_axes, env_df)
#
#   permanova_partial <- vegan::adonis2(
#     jaccard ~ env_label_good + PCoA1 + PCoA2 + PCoA3 + PCoA4 + PCoA5,
#     data = covariates,
#     permutations = 999,
#     by = "margin" # tests each term after all others — use for partial tests
#   ) |>
#     rownames_to_column("variable") |>
#     as.data.frame() |>
#     mutate(type = cat)
#
#   permanova_partial
# }

varpart_gene_phylogeny <- function(jaccard, tree, meta, cat) {
  genome_ids <- rownames(jaccard)

  phylo_dist <- cophenetic.phylo(tree@phylo)
  phylo_dist <- as.dist(phylo_dist[genome_ids, genome_ids])

  phylo_pcoa <- cmdscale(phylo_dist, k = 10, eig = TRUE)

  phylo_axes <- as.data.frame(phylo_pcoa$points[, 1:5])
  colnames(phylo_axes) <- paste0("PCoA", 1:5)

  env_df <- meta |>
    select(Genome, env_label_good) |>
    as.data.frame() |>
    arrange(match(Genome, genome_ids)) |>
    column_to_rownames("Genome")

  covariates <- cbind(phylo_axes, env_df)

  jaccard_mat <- as.matrix(jaccard)

  vp <- varpart(jaccard_mat, env_df["env_label_good"], phylo_axes)

  df <- vp$part$indfract |>
    as.data.frame() |>
    rownames_to_column("variables") |>
    mutate(
      variables = case_when(
        variables == "[a] = X1|X2" ~ "Environment",
        variables == "[b] = X2|X1" ~ "Phylogeny",
        variables == "[c]" ~ "Shared fraction",
        variables == "[d] = Residuals" ~ "Residuals"
      )
    ) |>
    select(-c("R.squared", "Testable")) |>
    mutate(type = cat)

  rda_env <- dbrda(
    jaccard ~ env_label_good +
      Condition(PCoA1 + PCoA2 + PCoA3 + PCoA4 + PCoA5),
    data = covariates
  )

  anova_rda_env <- as.data.frame(anova(rda_env)$`Pr(>F)`[1])

  rda_phy <- dbrda(
    jaccard ~ Condition(env_label_good) +
      PCoA1 +
      PCoA2 +
      PCoA3 +
      PCoA4 +
      PCoA5,
    data = covariates
  )

  anova_rda_phy <- as.data.frame(anova(rda_phy)$`Pr(>F)`[1])

  anova_rda_env_df <- anova_rda_env |>
    mutate(variables = "Environment") |>
    rename(p_value = "anova(rda_env)$`Pr(>F)`[1]")

  anova_rda_phy_df <- anova_rda_phy |>
    mutate(variables = "Phylogeny") |>
    rename(p_value = "anova(rda_phy)$`Pr(>F)`[1]")

  anova_rda_df <- rbind(anova_rda_env_df, anova_rda_phy_df)

  df_varpart <- df |>
    left_join(anova_rda_df, join_by(variables))

  df_varpart
}

adonis_gene_completeness <- function(jaccard, quality, meta, cat) {
  genome_ids <- rownames(jaccard)

  env_df <- meta |>
    select(Genome, Environment, Completeness) |>
    as.data.frame() |>
    arrange(match(Genome, genome_ids)) |>
    column_to_rownames("Genome")

  permanova <- adonis2(
    jaccard ~ Environment +
      Completeness,
    data = env_df,
    permutations = 999,
    by = "margin"
  ) |>
    rownames_to_column("variable") |>
    as.data.frame() |>
    mutate(type = cat)

  permanova
}

perturb_states <- function(states, p = 0.1) {
  new_states <- states
  idx <- sample(seq_along(states), size = round(p * length(states)))

  new_states[idx] <- sample(unique(states), length(idx), replace = TRUE)

  return(new_states)
}

run_ASR_tip <- function(tree, states) {
  fit <- ace(states, tree, type = "discrete", model = "ER")
  apply(fit$lik.anc, 1, which.max)
}

run_ace_tree <- function(tree, meta) {
  tree_real <- tree@phylo

  habitat <- meta |>
    select(Genome, env_label_good) |>
    arrange(match(Genome, tree_real$tip.label))

  ace <- ace(habitat$env_label_good, tree_real, type = "discrete", model = "ER")

  ace
}

sensitivity_tips <- function(ace, tree, meta, rate) {
  tree_real <- tree@phylo

  habitat <- meta |>
    select(Genome, env_label_good) |>
    arrange(match(Genome, tree_real$tip.label))

  original <- apply(ace[["lik.anc"]], 1, which.max)

  perturb_results <- replicate(100, {
    perturbed <- perturb_states(habitat$env_label_good, rate)
    run_ASR_tip(tree_real, perturbed)
  })

  agreement <- rowMeans(perturb_results == original)

  df <- data.frame(
    Rate = rate,
    Mean_agreement = round(mean(agreement, na.rm = TRUE), 3),
    SD_agreement = round(sd(agreement, na.rm = TRUE), 3),
    Min_agreement = round(min(agreement, na.rm = TRUE), 3)
  )

  df
}

get_ace_node_prob <- function(ace, meta, nodes, nodes_imp) {
  node_probs <- ace[["lik.anc"]]
  colnames(node_probs) <- colnames(ace[["lik.anc"]])
  max_prob <- apply(node_probs, 1, max)
  dominant <- colnames(node_probs)[apply(node_probs, 1, which.max)]

  conf_summary <- data.frame(
    Top_state = dominant,
    Max_prob = round(max_prob, 3),
    Entropy = round(
      apply(node_probs, 1, function(p) {
        p <- p[p > 0]
        -sum(p * log(p))
      }),
      3
    )
  ) |>
    rownames_to_column("label") |>
    left_join(nodes, join_by("label")) |>
    left_join(nodes_imp, join_by("node")) |>
    mutate(reason = if_else(is.na(reason), "Other", reason))

  conf_summary
}

get_max_prob <- function(prob, node_tree, nodes_imp) {
  df <- prob |>
    pivot_longer(cols = !node) |>
    group_by(node) |>
    slice(which.max(value)) |>
    rename(label = node) |>
    inner_join(node_tree, join_by("label")) |>
    left_join(nodes_imp, join_by("node")) |>
    mutate(reason = if_else(is.na(reason), "Other", reason))

  df
}

stochastic_mapping_acr <- function(tree, meta, past) {
  tree_real <- tree@phylo

  habitat <- meta |>
    select(Genome, env_label_good) |>
    arrange(match(Genome, tree_real$tip.label))

  tips <- habitat$env_label_good
  names(tips) <- tree_real$tip.label

  nsim <- 200
  smap <- make.simmap(
    tree_real,
    tips,
    model = "ER",
    nsim = nsim,
    message = FALSE
  )

  sm_sum <- summary(smap, plot = FALSE)

  ## Posterior probability of each state at each node
  simmap_pp <- sm_sum$ace
  colnames(simmap_pp) <- colnames(sm_sum$ace)

  simmap_pp_clean <- simmap_pp |>
    as.data.frame() |>
    rownames_to_column("type") |>
    filter(!str_detect(type, "_")) |>
    column_to_rownames("type")

  simmap_dominant <- colnames(simmap_pp_clean)[
    apply(simmap_pp_clean, 1, which.max)
  ]
  simmap_max_prob <- apply(simmap_pp_clean, 1, max)

  past_node <- past |>
    ungroup() |>
    select(node, name) |>
    arrange(node)

  agree <- mean(past_node$name == simmap_dominant)
  # cat(sprintf("ML vs SIMMAP state agreement: %.1f%%\n", 100 * agree))

  simmap_df <- data.frame(
    ML_state = past_node$name,
    SIMMAP_state = simmap_dominant,
    SIMMAP_pp = round(simmap_max_prob, 3)
  )
}

get_micro_cat <- function(micro) {
  list <- unique(micro$name)
}

run_caper_comp <- function(micro, meta, tree) {
  tree_clean <- tree@phylo

  meta_sel <- meta |>
    select(Genome, Completeness, Size)

  dat <- micro |>
    select(Genome, env_label_good, name, value) |>
    left_join(meta_sel, join_by(Genome)) |>
    mutate(log_genome_size = log(Size)) |>
    mutate(compl_z = as.numeric(scale(Completeness)))

  dat_wide <- dat |>
    pivot_wider(
      names_from = name,
      values_from = value,
      values_fill = 0
    ) |>
    as.data.frame()

  contrasts(dat_wide$env_label_good) <- contr.sum(nlevels(
    dat_wide$env_label_good
  ))

  comp_dat <- comparative.data(
    phy = tree_clean,
    data = dat_wide,
    names.col = "Genome",
    vcv = TRUE,
    warn.dropped = TRUE
  )

  comp_dat
}

fit_one_caper <- function(comp_dat, pw) {
  contrasts(comp_dat$data$env_label_good) <- contr.sum(nlevels(
    comp_dat$data$env_label_good
  ))

  fmla_full <- as.formula(
    paste(pw, "~ env_label_good + log_genome_size + compl_z")
  )
  fmla_null <- as.formula(
    paste(pw, "~ log_genome_size + compl_z")
  )

  full <- tryCatch(
    pgls(fmla_full, data = comp_dat, lambda = "ML"),
    error = function(e) NULL
  )
  null <- tryCatch(
    pgls(fmla_null, data = comp_dat, lambda = "ML"),
    error = function(e) NULL
  )

  if (is.null(full) || is.null(null)) {
    return(data.frame(
      pathway = pw,
      converged = FALSE,
      lambda = NA,
      LRT_stat = NA,
      df = NA,
      p_value = NA,
      stringsAsFactors = FALSE
    ))
  }

  ll_full <- logLik(full)
  ll_null <- logLik(null)
  lrt_stat <- 2 * (as.numeric(ll_full) - as.numeric(ll_null))
  lrt_df <- attr(ll_full, "df") - attr(ll_null, "df")
  p <- pchisq(lrt_stat, df = lrt_df, lower.tail = FALSE)

  data.frame(
    pathway = pw,
    converged = TRUE,
    lambda = round(full$param["lambda"], 4),
    LRT_stat = round(lrt_stat, 4),
    df = lrt_df,
    p_value = p,
    stringsAsFactors = FALSE
  )
}

clean_caper <- function(caper) {
  caper$signal_class <- with(
    caper,
    case_when(
      !converged ~ "failed",
      is.na(p_value) ~ "failed",
      p_value < 0.05 & lambda < 0.3 ~ "habitat_driven",
      p_value < 0.05 & lambda >= 0.7 ~ "confounded",
      p_value < 0.05 ~ "habitat_moderate_signal",
      p_value >= 0.05 & lambda >= 0.7 ~ "phylogenetically_conserved",
      p_value >= 0.05 & lambda < 0.3 ~ "neither",
      TRUE ~ "intermediate"
    )
  )

  caper <- caper[order(caper$p_value), ]
  caper
}

run_enrichment <- function(dat, count_col, label, n_boot = 5000) {
  dat$is_transition <- as.logical(dat$is_transition)

  x1 <- dat[[count_col]][dat$is_transition]
  x0 <- dat[[count_col]][!dat$is_transition]

  n_trans <- length(x1)
  n_back <- length(x0)

  ranks <- rank(dat[[count_col]], ties.method = "average")
  obs_stat <- mean(ranks[dat$is_transition])

  perm_p <- wilcox.test(
    x1,
    x0,
    alternative = "greater"
  )$p.value

  effect <- rank_biserial(x1, x0)$r_rank_biserial

  boot_eff <- replicate(n_boot, {
    s1 <- sample(x1, replace = TRUE)
    s0 <- sample(x0, replace = TRUE)
    rank_biserial(s1, s0)$r_rank_biserial
  })

  ci <- quantile(boot_eff, c(0.025, 0.975), na.rm = TRUE)

  list(
    label = label,
    perm_p = perm_p,
    effect = effect,
    ci_low = ci[1],
    ci_high = ci[2],
    obs_stat = obs_stat,
    n_trans = n_trans,
    n_back = n_back
  )
}

compare_nodes_transfers <- function(dtl, nodes, node_t, genomad) {
  events_count <- dtl |>
    select(species_label, duplications, losses, transfers, origination) |>
    right_join(node_t, join_by(species_label == label)) |>
    mutate(
      is_transition = node %in% nodes$node
    ) |>
    filter(!(node == 283))

  mge_type <- genomad |>
    filter(!is.na(species_label), !is.na(transfers), transfers > 0) |>
    group_by(type, species_label) |>
    summarise(
      mge_count = n(),
      mge_transfer_sum = sum(transfers, na.rm = TRUE),
      .groups = "drop"
    ) |>
    na.omit()

  results_list <- list()

  results_list[["HGT"]] <- run_enrichment(
    events_count,
    "transfers",
    "HGT - Alerax"
  )

  for (elem in c("plasmid", "virus", "all")) {
    if (elem == "all") {
      dat_mge <- events_count |>
        left_join(mge_type, join_by(species_label)) |>
        select(-type) |>
        group_by(species_label) |>
        summarise(
          across(where(is.numeric), sum, na.rm = TRUE),
          is_transition = any(is_transition),
          .groups = "drop"
        )
    } else {
      dat_mge <- events_count |>
        left_join(mge_type, join_by(species_label)) |>
        filter(type == elem)
    }

    label <- paste0("MGE — ", elem)
    results_list[[label]] <- run_enrichment(dat_mge, "mge_count", label)
  }

  summary_table <- map_dfr(results_list, function(res) {
    tibble(
      Test = res$label,
      `Permutation p` = signif(res$perm_p, 3),
      `Effect size (r)` = round(res$effect, 3),
      `CI low` = round(res$ci_low, 3),
      `CI high` = round(res$ci_high, 3),
      `Nodes trans/bg` = paste0(res$n_trans, " / ", res$n_back)
    )
  })
}
