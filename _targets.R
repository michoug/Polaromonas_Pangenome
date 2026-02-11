library(targets)
library(tarchetypes)
library(crew)
library(here)

set.seed(123)

tar_option_set(
  controller = crew_controller_local(
    workers = 4,
    seconds_idle = 5
  ),
  memory = "transient",
  garbage_collection = TRUE,
  # error = "continue",
  # cue = tar_cue(mode = "never")
)

source("R/packages.R")
source("R/functions.R")
source("R/plots.R")

# source("R/diverse.R")

set.seed(123)

tar_plan(
  tar_file_read(
    genome_stats,
    "rawData/genome_stats.tsv",
    read_tsv(!!.x, show_col_types = FALSE)
  ),

  tar_file_read(
    checkm2_data,
    "rawData/genome_quality_report.tsv",
    read_tsv(!!.x, show_col_types = FALSE)
  ),

  tar_file_read(
    genome_metadata,
    "rawData/genome_metadata.tsv",
    read_tsv(!!.x, show_col_types = FALSE)
  ),

  tar_file_read(
    genome_temp_grow,
    "rawData/genome_opT_growthRate.tsv",
    read_tsv(!!.x, show_col_types = FALSE)
  ),

  tar_file_read(
    gene_pres_abs,
    "rawData/gene_presence_absence.tsv.gz",
    read_tsv(!!.x, show_col_types = FALSE)
  ),

  tar_file_read(
    partition_persistent,
    "rawData/ppanggolin_persistent.txt.gz",
    read_tsv(!!.x, show_col_types = FALSE)
  ),

  tar_file_read(
    partition_shell,
    "rawData/ppanggolin_shell.txt.gz",
    read_tsv(!!.x, show_col_types = FALSE)
  ),

  tar_file_read(
    partition_cloud,
    "rawData/ppanggolin_cloud.txt.gz",
    read_tsv(!!.x, show_col_types = FALSE)
  ),

  tar_file_read(
    gene_families,
    "rawData/ppanggolin_gene_families.tsv.gz",
    read_tsv(!!.x, show_col_types = FALSE)
  ),

  tar_file_read(
    rooted_tree,
    "rawData/pastML_rooted_tree.nwk",
    read.nhx(!!.x)
  ),

  tar_file_read(
    dtl_all,
    "rawData/alerax_DTOL.csv",
    read_csv(!!.x, show_col_types = FALSE)
  ),

  tar_file_read(
    microtrait_data,
    "rawData/genome_microtrait_results.RDS",
    readRDS(!!.x)
  ),

  tar_file_read(rule2trait, "rawData/DB/rule2trait.RDS", readRDS(!!.x)),

  tar_file_read(allrules, "rawData/DB/rules.RDS", readRDS(!!.x)),

  tar_file_read(
    kegg_path,
    "rawData/DB/kegg_path.tsv.gz",
    read_tsv(!!.x, show_col_types = FALSE)
  ),

  tar_file_read(
    kegg_mod,
    "rawData/DB/kegg_modules.tsv.gz",
    read_tsv(!!.x, show_col_types = FALSE)
  ),

  tar_file_read(
    kegg_mod_prok,
    "rawData/DB/kegg_module_prok_5.txt",
    read_tsv(!!.x, show_col_types = FALSE)
  ),

  tar_file_read(
    eggnog_dat,
    "rawData/genome_eggnog.tsv.gz",
    read_tsv(!!.x, show_col_types = FALSE)
  ),

  tar_file_read(
    nodes_interest,
    "rawData/nodes_of_interest.txt",
    read_tsv(!!.x, show_col_types = FALSE)
  ),

  tar_file_read(
    rarefaction,
    "rawData/ppanggolin_rarefaction.csv.gz",
    read_csv(!!.x, show_col_types = FALSE)
  ),

  tar_target(
    events_all,
    combine_events_df("rawData/alerax_events_all.zip", "temp")
  ),

  tar_target(
    rds_families,
    combine_rds_df("rawData/gene_families_microtrait.zip", "temp2")
  ),

  tar_target(atlas_type, combine_atlas_df("rawData/atlas.zip", "temp3")),

  tar_target(
    genomad_virus,
    get_genomad_dat("rawData/genomad_virus.zip", "temp4", "cdcldddlddc")
  ),

  tar_target(
    genomad_plasmid,
    get_genomad_dat("rawData/genomad_plasmid.zip", "temp5", "cdcdddlddll")
  ),

  tar_target(
    proteins_contigs,
    get_contig_prot("rawData/polaromonas_contigs_proteins.zip", "temp6")
  ),

  tar_file_read(
    atlas_coord,
    "rawData/microbeAtlasCoordinates.txt",
    read_tsv(!!.x, show_col_types = FALSE)
  ),

  tar_file_read(
    atlas_samples,
    "rawData/microbeAtlasSamples.txt.gz",
    read_tsv(!!.x, show_col_types = FALSE)
  ),

  tar_file_read(
    collapse_info,
    "rawData/collapse_node.txt",
    read_tsv(!!.x, show_col_types = FALSE)
  ),

  tar_target(
    genome_metadata_comb,
    combine_metadata(genome_metadata, genome_stats, checkm2_data)
  ),

  tar_target(
    atlas_map,
    atlas_DF_map(atlas_type, atlas_coord, atlas_samples)
  ),

  tar_target(genome_metadata_clean, clean_env_metadata(genome_metadata)),

  tar_target(colors_samples, getColors()),

  tar_target(colors_countries, getColorsCountry(genome_metadata_clean)),

  tar_target(colors_DTL, getcolorsDTL()),

  tar_target(
    genome_stat_cat,
    combine_stats_checkm_temp(genome_stats, checkm2_data, genome_temp_grow)
  ),

  tar_target(
    genome_statToPlot,
    combine_with_groups(genome_stat_cat, genome_metadata_clean)
  ),

  tar_target(
    nmds_persistent_dat,
    prepare_genes_for_nmds(gene_pres_abs, partition_persistent)
  ),

  tar_target(nmds_persistent, create_nmds(nmds_persistent_dat)),

  tar_target(
    anosim_persistent,
    anosim_nmds(nmds_persistent_dat, genome_metadata_clean)
  ),

  tar_target(
    nmds_shell_dat,
    prepare_genes_for_nmds(gene_pres_abs, partition_shell)
  ),

  tar_target(nmds_shell, create_nmds(nmds_shell_dat)),

  tar_target(anosim_shell, anosim_nmds(nmds_shell_dat, genome_metadata_clean)),

  tar_target(
    nmds_cloud_dat,
    prepare_genes_for_nmds(gene_pres_abs, partition_cloud)
  ),

  tar_target(nmds_cloud, create_nmds(nmds_cloud_dat)),

  tar_target(anosim_cloud, anosim_nmds(nmds_cloud_dat, genome_metadata_clean)),

  tar_target(
    data_env_stats,
    clean_env_stats(genome_metadata_clean, genome_stat_cat)
  ),

  tar_target(
    microtrait_to_plot,
    prepare_microtrait_data(microtrait_data, data_env_stats)
  ),

  tar_target(
    microtrait_enrich_sign,
    stats_microtrait_enrich(microtrait_to_plot)
  ),

  tar_target(eggnog_ko, clean_egg(eggnog_dat)),

  tar_target(families_kos, get_families_Kos(gene_families, eggnog_ko)),

  tar_target(events_clean, clean_events(events_all)),

  tar_target(node_tree, get_node_tree(rooted_tree)),

  tar_target(events_ko, events_kos(events_clean, families_kos, node_tree)),

  tar_target(
    DTL_nodes_int,
    getDTLNodes(events_clean, node_tree, nodes_interest)
  ),

  tar_target(rule2trait_clean, clean_rules(rule2trait, allrules)),

  tar_target(
    events_micro,
    get_micro_events(events_clean, rds_families, rule2trait_clean, node_tree)
  ),

  tar_target(kegg_mod_clean, clean_kegg_mod(kegg_mod, kegg_mod_prok)),

  tar_target(kegg_path_clean, clean_kegg_path(kegg_path)),

  tar_target(
    ev_kegg_ancester,
    getFunctionAncesterKO(events_ko, kegg_path_clean, kegg_mod_clean),
  ),

  tar_target(
    ev_kegg_ancester_sel,
    selectFunctionAncester(ev_kegg_ancester, nodes_interest),
    pattern = map(nodes_interest)
  ),

  tar_target(ev_micro_ancester, getFunctionAncesterMicro(events_micro)),

  tar_target(
    ev_micro_ancester_sel,
    selectFunctionAncester(ev_micro_ancester, nodes_interest),
    pattern = map(nodes_interest)
  ),

  tar_target(
    figure_map_atlas,
    plot_map_atlas(atlas_map)
  ),

  tar_target(
    figure_map_atlas_pdf,
    ggsave(
      "Figures/Figure_S1_map_atlas.pdf",
      figure_map_atlas,
      width = 12,
      height = 7.5,
      create.dir = T
    ),
    format = "file"
  ),

  tar_target(
    figure_map_atlas_png,
    ggsave(
      "Figures/Figure_S1_map_atlas.png",
      figure_map_atlas,
      width = 12,
      height = 7.5,
      dpi = 300,
      create.dir = T
    ),
    format = "file"
  ),

  tar_target(
    map_data_world,
    prepare_map_data(genome_metadata_clean, colors_samples, "world")
  ),

  tar_target(
    map_data_eur,
    prepare_map_data(genome_metadata_clean, colors_samples, "europe")
  ),

  tar_target(
    figure_map_world,
    plot_map(map_data_world, colors_samples, "world")
  ),

  tar_target(
    figure_map_eur,
    plot_map(map_data_eur, colors_samples, "europe")
  ),

  tar_target(table_stats, draw_table(genome_statToPlot)),

  tar_target(
    table_meta_tsv,
    write_tsv(genome_metadata_comb, "Figures/Table_S1_meta.tsv")
  ),

  tar_target(
    table_stats_docx,
    gtsave(table_stats, "Figures/Table_S2_stats.docx"),
    format = "file"
  ),

  tar_target(
    table_stats_pdf,
    gtsave(table_stats, "Figures/Table_S2_stats.pdf"),
    format = "file",
    packages = c("gt", "webshot2")
  ),

  tar_target(
    figure_genome_comparison,
    comparing_genomes_params(genome_statToPlot, colors_samples)
  ),

  tar_target(
    figure_genome_comparison_pdf,
    ggsave(
      "Figures/Figure_1b_comparison.pdf",
      figure_genome_comparison,
      width = 12,
      height = 8,
      create.dir = T
    ),
  ),

  tar_target(
    figure_genome_comparison_png,
    ggsave(
      "Figures/Figure_1b_comparison.png",
      figure_genome_comparison,
      width = 12,
      height = 8,
      dpi = 300,
      create.dir = T
    ),
    format = "file"
  ),

  tar_target(
    figure_1_map_comp,
    p <- (figure_map_world | figure_map_eur) /
      figure_genome_comparison +
      plot_layout(guides = "collect", heights = c(1, 0.5)) +
      plot_annotation(tag_levels = c("A", "", "B"))
  ),

  tar_target(
    figure_1_map_comp_pdf,
    ggsave(
      "Figures/Figure_1.pdf",
      figure_1_map_comp,
      width = 12,
      height = 8,
      create.dir = T
    )
  ),

  tar_target(
    figure_1_map_comp_png,
    ggsave(
      "Figures/Figure_1.png",
      figure_1_map_comp,
      width = 12,
      height = 8,
      create.dir = T,
      dpi = 300
    ),
  ),

  tar_target(
    figure_rarefaction,
    rarefaction_curve(rarefaction)
  ),

  tar_target(
    figure_rarefaction_pdf,
    ggsave(
      "Figures/Figure_S2_rarefaction_curve.pdf",
      figure_rarefaction,
      width = 10,
      height = 6,
      create.dir = T
    )
  ),

  tar_target(
    figure_rarefaction_png,
    ggsave(
      "Figures/Figure_S2_rarefaction_curve.png",
      figure_rarefaction,
      width = 10,
      height = 6,
      dpi = 300,
      create.dir = T
    )
  ),

  tar_target(
    figure_S3_prevalence_ppan,
    plot_prevalence_ppan(
      gene_pres_abs,
      partition_persistent,
      partition_shell,
      partition_cloud
    )
  ),

  tar_target(
    figure_S3_prevalence_ppan_pdf,
    ggsave(
      "Figures/Figure_S3_prevalence_ppan.pdf",
      figure_S3_prevalence_ppan,
      width = 8,
      height = 6,
      create.dir = T
    )
  ),

  tar_target(
    figure_S3_prevalence_ppan_png,
    ggsave(
      "Figures/Figure_S3_prevalence_ppan.png",
      figure_S3_prevalence_ppan,
      width = 8,
      height = 6,
      create.dir = T
    )
  ),

  tar_target(
    figure_nmds_persistent,
    plot_nmds(nmds_persistent, data_env_stats, colors_samples)
  ),

  tar_target(
    figure_nmds_persistent_pdf,
    ggsave(
      "Figures/Figure_2_NMDS_persistent.pdf",
      figure_nmds_persistent,
      width = 8,
      height = 8,
      create.dir = T
    )
  ),

  tar_target(
    figure_nmds_persistent_png,
    ggsave(
      "Figures/Figure_2_NMDS_persistent.png",
      figure_nmds_persistent,
      width = 8,
      height = 8,
      dpi = 300,
      create.dir = T
    )
  ),

  tar_target(
    figure_nmds_shell,
    plot_nmds(nmds_shell, data_env_stats, colors_samples)
  ),

  tar_target(
    figure_nmds_cloud,
    plot_nmds(nmds_cloud, data_env_stats, colors_samples)
  ),

  tar_target(
    figure_nmds_cloud_shell,
    wrap_plots(
      figure_nmds_shell,
      figure_nmds_cloud,
      guides = "collect",
      nrow = 1
    ) +
      plot_annotation(tag_levels = 'A') &
      theme(legend.position = 'bottom')
  ),

  tar_target(
    figure_nmds_cloud_shell_pdf,
    ggsave(
      "Figures/Figure_S4_NMDS_shell_cloud.pdf",
      figure_nmds_cloud_shell,
      width = 15,
      height = 10,
      create.dir = T
    )
  ),

  tar_target(
    figure_nmds_cloud_shell_png,
    ggsave(
      "Figures/Figure_S4_NMDS_shell_cloud.png",
      figure_nmds_cloud_shell,
      width = 15,
      height = 10,
      dpi = 300,
      create.dir = T
    )
  ),

  tar_target(
    figure_heatmap_microtrait,
    plot_heatmap_microtrait(microtrait_to_plot, microtrait_enrich_sign)
  ),

  # tar_target(
  #   figure_heatmap_microtrait_pdf,
  #   ggsave(
  #     "Figures/Figure_X_heatmap_microtrait.pdf",
  #     figure_heatmap_microtrait,
  #     width = 12,
  #     height = 8,
  #     create.dir = T
  #   )
  # ),

  tar_target(
    collapse_color,
    collapse_node_color(collapse_info, colors_samples)
  ),

  tar_target(
    figure_phylogenetic_tree,
    phylogenetic_tree(
      rooted_tree,
      nodes_interest,
      colors_samples,
      collapse_color
    )
  ),

  tar_target(
    figure_phylogenetic_tree_pdf,
    ggsave(
      "Figures/Figure_3_phylogenetic_tree.pdf",
      figure_phylogenetic_tree,
      width = 15,
      height = 10,
      create.dir = T
    )
  ),

  tar_target(
    figure_phylogenetic_tree_png,
    ggsave(
      "Figures/Figure_3_phylogenetic_tree.png",
      figure_phylogenetic_tree,
      width = 15,
      height = 10,
      dpi = 300,
      create.dir = T
    )
  ),

  tar_target(
    figure_DTL_nodes,
    plot_info_dtl(DTL_nodes_int, colors_DTL),
  ),

  tar_target(
    figure_DTL_nodes_pdf,
    ggsave(
      "Figures/Figure_S5_DTL_nodes.pdf",
      figure_DTL_nodes,
      width = 12,
      height = 8,
      create.dir = T
    )
  ),

  tar_target(
    figure_DTL_nodes_png,
    ggsave(
      "Figures/Figure_S5_DTL_nodes.png",
      figure_DTL_nodes,
      width = 12,
      height = 8,
      dpi = 300,
      create.dir = T
    )
  ),

  tar_target(
    figure_heatmap_DTL,
    heatmap_micro_dtl(ev_micro_ancester_sel)
  ),

  tar_target(
    figure_heatmap_DTL_pdf,
    ggsave(
      "Figures/Figure_4_heatmap_DTL.pdf",
      figure_heatmap_DTL,
      width = 15,
      height = 12,
      create.dir = T
    )
  ),

  tar_target(
    figure_heatmap_DTL_png,
    ggsave(
      "Figures/Figure_4_heatmap_DTL.png",
      figure_heatmap_DTL,
      width = 15,
      height = 12,
      dpi = 300,
      create.dir = T
    )
  ),

  tar_target(
    genomad_dtl,
    get_genomad_dtl(
      proteins_contigs,
      gene_families,
      events_clean,
      genomad_virus,
      genomad_plasmid
    )
  ),

  tar_target(
    genomad_numb,
    genomad_bar(genome_metadata_clean, genomad_virus, genomad_plasmid)
  ),

  tar_target(
    figure_genomad_numb,
    plot_genomad_perc(genomad_numb)
  ),

  # tar_target(
  #   figure_genomad_numb_pdf,
  #   ggsave(
  #     "Figures/Figure_S6a_Virus_Plasmid.pdf",
  #     figure_genomad_numb,
  #     width = 12,
  #     height = 8,
  #     create.dir = T
  #   )
  # ),
  #
  # tar_target(
  #   figure_genomad_numb_png,
  #   ggsave(
  #     "Figures/Figure_S6a_Virus_Plasmid.png",
  #     figure_genomad_numb,
  #     width = 12,
  #     height = 8,
  #     dpi = 300,
  #     create.dir = T
  #   )
  # ),

  tar_target(
    genomad_dtl_reason,
    genomad_dtl_filter(genomad_dtl, node_tree, nodes_interest)
  ),

  tar_target(
    figure_genomad_dtl,
    plot_genomad_dtl(genomad_dtl_reason, colors_DTL)
  ),

  # tar_target(
  #   figure_genomad_dtl_pdf,
  #   ggsave(
  #     "Figures/Figure_S6b_Virus_Plasmid_DTL.pdf",
  #     figure_genomad_dtl,
  #     width = 12,
  #     height = 8,
  #     create.dir = T
  #   )
  # ),
  #
  # tar_target(
  #   figure_genomad_dtl_png,
  #   ggsave(
  #     "Figures/Figure_S6b_Virus_Plasmid_DTL.png",
  #     figure_genomad_dtl,
  #     width = 12,
  #     height = 8,
  #     dpi = 300,
  #     create.dir = T
  #   )
  # ),

  tar_target(
    figure_genomad,
    wrap_plots(
      figure_genomad_numb,
      figure_genomad_dtl,
      ncol = 1
    ) +
      plot_annotation(tag_levels = 'A') &
      theme(legend.position = 'bottom')
  ),

  tar_target(
    figure_genomad_pdf,
    ggsave(
      "Figures/Figure_S6_Virus_Plasmid.pdf",
      figure_genomad,
      width = 15,
      height = 11,
      create.dir = T
    )
  ),

  tar_target(
    figure_genomad_png,
    ggsave(
      "Figures/Figure_S6_Virus_Plasmid.png",
      figure_genomad,
      width = 15,
      height = 11,
      dpi = 300,
      create.dir = T
    )
  ),
)

# targets::tar_make(reporter = "balanced")
