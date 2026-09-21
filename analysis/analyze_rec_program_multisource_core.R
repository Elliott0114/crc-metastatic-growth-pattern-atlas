#!/usr/bin/env Rscript

# Build an interpretable cross-dataset recurrence map for the fixed REC programme.

suppressPackageStartupMessages({
  library(data.table)
})

options(stringsAsFactors = FALSE)

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
root <- if (length(file_arg)) {
  normalizePath(file.path(dirname(sub("^--file=", "", file_arg[[1]])), ".."), mustWork = TRUE)
} else {
  normalizePath(".", mustWork = TRUE)
}

out_dir <- file.path(root, "analysis_results", "rec_program_multisource_core")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_result <- function(...) {
  fread(file.path(root, ...))
}

write_tsv <- function(x, filename) {
  fwrite(as.data.table(x), file.path(out_dir, filename), sep = "\t", na = "NA")
}

program <- read_result(
  "analysis_results", "ogden_anchor_program", "selected_anchor_program_top50.tsv"
)[, .(gene = toupper(gene), program_rank, discovery_rec_vs_other = logFC)]

bulk <- read_result(
  "analysis_results", "rec_program_gse151165_bulk_projection", "gene_level_effects.tsv"
)[, .(gene = toupper(gene), bulk_rHGP_vs_dHGP = rhgp_minus_dhgp)]

interface <- read_result(
  "analysis_results", "rec_program_latacz_interface_summary", "gene_effects.tsv"
)[gene_set == "frozen_rec_top50",
  .(gene = toupper(gene), interface_rHGP_vs_dHGP = uncorrected_log2_r_over_d)]

epithelial <- read_result(
  "analysis_results", "rec_program_cross_modal_concordance", "gene_level_concordance.tsv"
)[, .(
  gene = toupper(gene),
  epithelial_rHGP_vs_dHGP = scrna_rhgp_minus_dhgp,
  spatial_rHGP_vs_dHGP = spatial_rhgp_minus_dhgp
)]

outgrowth <- read_result(
  "analysis_results", "gse294385_rec_program_extension", "gene_level_paired_effects.tsv"
)[as.character(member_rec_equal_gene) %chin% c("True", "TRUE", "true"),
  .(gene = toupper(gene), outgrowth_macro_vs_micro = mean_macro_minus_micro)]

modules <- read_result(
  "analysis_results", "rec_program_two_axis_interpretation", "gene_two_axis_classification.tsv"
)[, .(
  gene = toupper(gene),
  member_ap1 = as.character(member_experimental_ap1_target_overlap) %chin% c("True", "TRUE", "true"),
  member_hypoxia = as.character(member_hypoxia_mp6_overlap) %chin% c("True", "TRUE", "true"),
  member_regenerative = as.character(member_regenerative_overlap) %chin% c("True", "TRUE", "true"),
  member_nfkb = as.character(member_nfkb_regulon_overlap) %chin% c("True", "TRUE", "true")
)]

gene_map <- Reduce(
  function(x, y) merge(x, y, by = "gene", all = FALSE),
  list(program, bulk, interface, epithelial, outgrowth, modules)
)

replication_columns <- c(
  "bulk_rHGP_vs_dHGP", "interface_rHGP_vs_dHGP",
  "epithelial_rHGP_vs_dHGP", "spatial_rHGP_vs_dHGP",
  "outgrowth_macro_vs_micro"
)
hgp_columns <- replication_columns[1:4]

gene_map[, hgp_positive_sources := rowSums(.SD > 0), .SDcols = hgp_columns]
gene_map[, all_positive_sources := rowSums(.SD > 0), .SDcols = replication_columns]
gene_map[, persistent_five_source_core := all_positive_sources == length(replication_columns)]
gene_map[, hgp_recurrence_class := fcase(
  hgp_positive_sources == 4L, "positive in all four HGP datasets",
  hgp_positive_sources == 3L, "positive in three HGP datasets",
  hgp_positive_sources == 2L, "split HGP direction",
  default = "positive in zero or one HGP dataset"
)]
gene_map[, outgrowth_direction := fifelse(
  outgrowth_macro_vs_micro > 0, "higher in macrometastasis", "not higher in macrometastasis"
)]
gene_map[, biological_membership := apply(.SD, 1L, function(row) {
  labels <- c("AP-1 target", "Hypoxia", "Regenerative", "NF-kB regulon")[as.logical(row)]
  if (!length(labels)) "Other REC programme" else paste(labels, collapse = "; ")
}), .SDcols = c("member_ap1", "member_hypoxia", "member_regenerative", "member_nfkb")]

effect_long <- melt(
  gene_map,
  id.vars = c(
    "gene", "program_rank", "persistent_five_source_core", "hgp_positive_sources",
    "all_positive_sources", "hgp_recurrence_class", "outgrowth_direction",
    "biological_membership"
  ),
  measure.vars = c("discovery_rec_vs_other", replication_columns),
  variable.name = "source",
  value.name = "log2_effect"
)
source_labels <- c(
  discovery_rec_vs_other = "REC discovery",
  bulk_rHGP_vs_dHGP = "Bulk HGP",
  interface_rHGP_vs_dHGP = "Interface HGP",
  epithelial_rHGP_vs_dHGP = "Epithelial HGP",
  spatial_rHGP_vs_dHGP = "FFPE spatial HGP",
  outgrowth_macro_vs_micro = "Macro vs micro"
)
effect_long[, source_label := source_labels[as.character(source)]]
effect_long[, display_effect := pmax(-2, pmin(2, log2_effect))]

source_summary <- rbindlist(lapply(replication_columns, function(column) {
  data.table(
    source = column,
    source_label = source_labels[[column]],
    n_genes = nrow(gene_map),
    positive_n = sum(gene_map[[column]] > 0),
    positive_fraction = mean(gene_map[[column]] > 0),
    median_log2_effect = median(gene_map[[column]]),
    mean_log2_effect = mean(gene_map[[column]])
  )
}))

core <- gene_map[persistent_five_source_core == TRUE][order(-outgrowth_macro_vs_micro)]

correlation_matrix <- cor(
  gene_map[, ..replication_columns],
  method = "spearman",
  use = "pairwise.complete.obs"
)
correlation_long <- as.data.table(as.table(correlation_matrix))
setnames(correlation_long, c("source_1", "source_2", "spearman_rho"))
correlation_long[, source_1_label := source_labels[as.character(source_1)]]
correlation_long[, source_2_label := source_labels[as.character(source_2)]]

write_tsv(gene_map[order(program_rank)], "gene_recurrence_map.tsv")
write_tsv(effect_long[order(program_rank, source)], "gene_effects_long.tsv")
write_tsv(source_summary, "source_summary.tsv")
write_tsv(core, "persistent_core.tsv")
write_tsv(correlation_long, "source_correlations.tsv")

run_info <- data.table(
  analysis_date = "2026-08-23",
  method = "transparent direction count across four HGP datasets plus paired outgrowth",
  n_common_genes = nrow(gene_map),
  n_persistent_core = nrow(core),
  persistent_core_rule = "positive effect in bulk HGP, interface HGP, epithelial HGP, FFPE spatial HGP, and paired macro-versus-micro data",
  r_version = R.version.string,
  data_table_version = as.character(packageVersion("data.table"))
)
write_tsv(run_info, "run_info.tsv")
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))

summary_lines <- c(
  "# Multi-source REC-program recurrence map",
  "",
  sprintf(
    "%d programme genes were measurable in all four HGP datasets and the paired outgrowth dataset.",
    nrow(gene_map)
  ),
  paste0(
    "Nine genes were positive in every replication context: ",
    paste(core$gene, collapse = ", "), "."
  ),
  "",
  paste0(
    "The recurrence rule is a transparent sign count. It does not treat datasets as exchangeable replicates, ",
    "fit a latent trajectory or assign causal order."
  )
)
writeLines(summary_lines, file.path(out_dir, "analysis_summary.md"), useBytes = TRUE)

manifest <- data.table(
  file = c(
    "gene_recurrence_map.tsv", "gene_effects_long.tsv", "source_summary.tsv",
    "persistent_core.tsv", "source_correlations.tsv", "run_info.tsv",
    "sessionInfo.txt", "analysis_summary.md"
  ),
  role = c(
    "gene-level recurrence classification", "long effect table for plotting", "source-level direction summary",
    "nine-gene persistent core", "descriptive source correlations", "method record",
    "R environment", "human-readable result summary"
  )
)
fwrite(manifest, file.path(out_dir, "_analysis_outputs.md"), sep = "\t")

