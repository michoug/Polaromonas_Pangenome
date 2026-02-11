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
    left_join(check, join_by("user_genome" == "Name")) %>%
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

  unzip(zipfile, exdir = folder, overwrite = T)

  lst <- dir_ls(folder) |>
    set_names(basename) |>
    map(read_tsv, show_col_types = F) |>
    list_rbind(names_to = "type")

  unlink(folder, recursive = T, force = T)
  lst
}

atlas_DF_map <- function(df, coord, samples) {
  df_clean <- df %>%
    clean_names() %>%
    mutate(samples = paste(run_id, number_sample_id, sep = ".")) %>%
    left_join(samples, join_by("samples")) %>%
    mutate(
      nbRead = abundance_ppm * RNA_reads / 1e6,
      percAbund = nbRead * 100 / RNA_reads
    ) %>%
    filter(percAbund > 0.1) %>%
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
    ) %>%
    mutate(env_good = str_to_sentence(env_good)) %>%
    group_by(number_sample_id, samples, env_good) %>%
    summarise(percAbund = sum(percAbund)) %>%
    ungroup() %>%
    left_join(coord, join_by(number_sample_id == sampleId)) %>%
    filter(!(parsedLon == "None")) %>%
    select(
      number_sample_id,
      samples,
      percAbund,
      env_good,
      parsedLat,
      parsedLon
    ) %>%
    filter(!(number_sample_id == "SRS7752820")) %>%
    mutate(parsedLon = as.numeric(parsedLon)) %>%
    mutate(parsedLat = as.numeric(parsedLat)) %>%
    filter(parsedLat < 90) %>%
    mutate(env_good = as.factor(env_good)) %>%
    mutate(env_good = factor(env_good, levels = sort(unique(env_good)))) %>%
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

  unzip(zipfile, exdir = folder, overwrite = T)

  lst <- dir_ls(folder) |>
    set_names(basename) |>
    map(
      read_tsv,
      show_col_types = F,
      progress = F,
      col_types = "cc",
      comment = "#",
      col_names = c("Contig", "Protein")
    ) |>
    list_rbind(names_to = "Genome") |>
    mutate(
      Genome = Genome |>
        str_replace_all(".tsv", "")
    )

  unlink(folder, recursive = T, force = T)
  lst
}

get_genomad_dtl <- function(prot, families, ev, virus, plasmid) {
  vir <- virus %>%
    mutate(type = "virus") %>%
    select(Genome, seq_name, type)

  plas <- plasmid %>%
    mutate(type = "plasmid") %>%
    select(Genome, seq_name, type)

  dat <- vir %>%
    rbind(plas) %>%
    select(-Genome) %>%
    left_join(
      prot,
      join_by("seq_name" == "Contig"),
      relationship = "many-to-many"
    ) %>%
    left_join(families, join_by("Protein" == "Gene")) %>%
    select(-c("localGene", "Frag", "Protein")) %>%
    distinct() %>%
    na.omit() %>%
    left_join(ev, join_by("Family"), relationship = "many-to-many")
}

prepare_genes_for_nmds <- function(genes, partition) {
  dat <- genes |>
    filter(Gene %in% partition$Family) |>
    column_to_rownames("Gene") |>
    select(where(function(x) sum(x) > 0)) |>
    as.matrix() |>
    t()
}

create_nmds <- function(dat) {
  nmds <- metaMDS(
    dat,
    k = 2,
    trace = 1,
    autotransform = F,
    distance = "jaccard",
    trymax = 100
  )

  nmds
}

anosim_nmds <- function(dat, map) {
  rownames(dat) <- gsub("-", "_", rownames(dat))

  map_clean <- map %>%
    arrange(match(Genome, rownames(dat)))

  dat_dist <- (dat) %>%
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
    mutate(value = if_else(value >= 1, 1, 0))

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

# get_families_Kos <- function(dat, egg){
#
#   # Apply majority rule to determine which Kos are in which families
#   dat_clean <- dat |>
#     left_join(egg, join_by(Gene == query))|>
#     filter(!is.na(KEGG_ko)) |>
#     group_by(Family)|>
#     mutate(nbGenes = n())|>
#     ungroup() |>
#     group_by(Family, KEGG_ko, nbGenes)|>
#     summarise(nbGenesPerKo = n())|>
#     mutate(perc = nbGenesPerKo * 100 / nbGenes)|>
#     ungroup()|>
#     group_by(Family) |>
#     mutate(nbKoperFam = n())|>
#     slice_max(perc, n = 1)|>
#     select(Family, KEGG_ko)
#
#   dat_clean
# }

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


get_node_tree <- function(tree) {
  p1 <- ggtree(tree)

  node_data <- p1$data |>
    filter(isTip == F) |>
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

# process_rds <- function(path) {
#   df <- readRDS(path)
#
#   df_clean <- df$genes_detected_table |>
#     mutate(name = gsub(".*/(.*).microtrait.rds", "\\1", path)) |>
#     select(name, gene_name, hmm_name) |>
#     as.data.frame()
# }

process_rds <- function(path) {
  df <- readRDS(path)

  df_clean <- df$rules_asserted |>
    clean_names() %>%
    mutate(name = gsub(".*/(.*).microtrait.rds", "\\1", path)) |>
    filter(microtrait_rule_asserted == T) %>%
    select(name, microtrait_rule_name, microtrait_rule_boolean) %>%
    mutate(
      parsed_components = map(microtrait_rule_boolean, parse_boolean_expression)
    ) %>%
    unnest(parsed_components) |>
    select(-microtrait_rule_boolean) %>%
    distinct() %>%
    as.data.frame()

  df_genes <- df$genes_detected_table %>%
    select(gene_name, hmm_name)

  df_dup <- df_clean %>%
    filter(microtrait_rule_name != parsed_components)

  list_dup <- intersect(df_dup$microtrait_rule_name, df_dup$parsed_components)

  df_clean_sel <- df_clean %>%
    select(-name)

  df_non_dup <- df_clean %>%
    full_join(
      df_clean_sel,
      join_by(parsed_components == microtrait_rule_name)
    ) %>%
    mutate(
      parsed_components = if_else(
        is.na(parsed_components.y),
        parsed_components,
        parsed_components.y
      )
    ) %>%
    select(-parsed_components.y) %>%
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

  unzip(zipfile, exdir = folder, overwrite = T)

  lst <- dir_ls(folder) |>
    map(process_rds) |>
    list_rbind()

  unlink(folder, recursive = T, force = T)
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
    ) %>%
    distinct() %>%
    clean_names()
}

subset_nodes <- function(ev, node, node_int) {
  dat <- ev %>%
    left_join(node, join_by(species_label == label)) %>%
    right_join(node_int, join_by(node))
}

get_micro_events <- function(ev, rds, rules, node) {
  rds_clean <- rules %>%
    mutate(microtrait_rule_name = as.character(microtrait_rule_name)) %>%
    right_join(
      rds,
      join_by(microtrait_rule_name),
      relationship = "many-to-many"
    ) %>%
    filter(!is.na(microtrait_trait_name3)) %>%
    select(-name) %>%
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
    ) %>%
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

  unzip(zipfile, exdir = folder, overwrite = T)

  lst <- dir_ls(folder) |>
    set_names(basename) |>
    map(
      read_tsv,
      show_col_types = F,
      progress = F,
      col_types = type
    ) |>
    list_rbind(names_to = "Genome") |>
    mutate(
      Genome = Genome |>
        str_replace_all("-", "_") |>
        str_remove("_(virus|plasmid)_summary\\.tsv$")
    )

  unlink(folder, recursive = T, force = T)
  lst
}

genomad_bar <- function(meta, virus, plasmid) {
  dat_virus <- virus %>%
    select(Genome) %>%
    mutate(type = "Virus") %>%
    distinct()

  dat_plasmid <- plasmid %>%
    select(Genome) %>%
    mutate(type = "Plasmid") %>%
    distinct()

  meta <- meta %>%
    group_by(env_label_good) %>%
    mutate(n_env = n()) %>%
    ungroup() %>%
    select(Genome, env_label_good, n_env)

  dat_all <- dat_virus %>%
    rbind(dat_plasmid) %>%
    left_join(meta, join_by("Genome")) %>%
    mutate(env_n = glue("{env_label_good} ({n_env})")) %>%
    group_by(type, env_label_good, env_n) %>%
    reframe(n = n() * 100 / n_env) %>%
    distinct()

  dat_all
}

genomad_dtl_filter <- function(genomad, tree, nodes) {
  genomad_dtl_filter <- genomad %>%
    left_join(tree, join_by("species_label" == "label")) %>%
    right_join(nodes, join_by("node")) %>%
    mutate(Genome = gsub("-", "_", Genome)) %>%
    group_by(type, reason) %>%
    mutate(seq_genome = paste(Genome, seq_name, sep = "_")) %>%
    mutate(allnumb = n_distinct(seq_genome)) %>%
    ungroup() %>%
    select(
      seq_genome,
      type,
      Family,
      duplications:reason,
      allnumb
    ) %>%
    distinct() %>%
    group_by(seq_genome, reason) %>%
    mutate(
      n = n(),
      dup_sum = sum(duplications) / n,
      trans_sum = sum(transfers) / n,
      orig_sum = sum(origination) / n,
      loss_sum = sum(losses) / n
    ) %>%
    ungroup() %>%
    select(seq_genome, type, reason, allnumb, n:loss_sum) %>%
    distinct() %>%
    rowwise() %>%
    mutate(max_row = max(c_across(dup_sum:loss_sum))) %>%
    filter(max_row > 0.2) %>%
    select(-c(max_row, n)) %>%
    pivot_longer(cols = !c(seq_genome, type, reason, allnumb)) %>%
    filter(value > 0) %>%
    group_by(seq_genome, type, reason, allnumb) %>%
    slice_max(order_by = value, n = 1, with_ties = F) %>%
    ungroup() %>%
    group_by(type, reason, allnumb, name) %>%
    reframe(perc = n() / allnumb) %>%
    distinct()
}
