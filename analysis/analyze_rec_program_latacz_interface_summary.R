#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
})

options(stringsAsFactors = FALSE)

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
root <- if (length(file_arg)) {
  normalizePath(file.path(dirname(sub("^--file=", "", file_arg[[1]])), ".."), mustWork = TRUE)
} else {
  normalizePath(".", mustWork = TRUE)
}

workbook <- file.path(
  root, "data_sources", "Latacz_2024_targeted_HGP_interface",
  "SupplTablesLataczforupload.xlsx"
)
program_path <- file.path(
  root, "analysis_results", "ogden_anchor_program", "selected_anchor_program_top50.tsv"
)
core_path <- file.path(
  root, "analysis_results", "rec_program_cross_modal_concordance",
  "top15_interpretable_core.tsv"
)
spec_path <- file.path(
  root, "metadata", "rec_program_latacz_interface_summary_spec_2026-08-23.md"
)
out_dir <- file.path(root, "analysis_results", "rec_program_latacz_interface_summary")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

blocks <- c("uncorrected", "treatment", "cms", "ts_ratio")
effect_columns <- paste0(blocks, "_log2_r_over_d")

write_tsv <- function(x, filename) {
  fwrite(as.data.table(x), file.path(out_dir, filename), sep = "\t", na = "NA")
}

raw <- as.data.table(read_excel(
  workbook,
  sheet = "SuppTab1",
  col_names = FALSE,
  skip = 3,
  .name_repair = "minimal"
))
raw <- raw[, seq_len(23L), with = FALSE]
columns <- c("ensembl", "gene", "entrez")
for (block in blocks) {
  columns <- c(
    columns,
    paste0(block, c("_expression", "_log2_r_over_d", "_se", "_p", "_fdr"))
  )
}
setnames(raw, columns)
raw[, gene := toupper(trimws(as.character(gene)))]
numeric_columns <- setdiff(columns, c("ensembl", "gene", "entrez"))
raw[, (numeric_columns) := lapply(
  .SD,
  function(values) suppressWarnings(as.numeric(values))
), .SDcols = numeric_columns]
raw <- unique(raw[!is.na(gene) & gene != ""], by = "gene")

program <- fread(program_path)
program[, gene := toupper(gene)]
core <- fread(core_path)
core[, gene := toupper(gene)]
gene_sets <- list(
  frozen_rec_top50 = program$gene,
  cross_modal_hgp_top15 = core$gene
)

summary_rows <- list()
gene_rows <- list()
for (gene_set in names(gene_sets)) {
  requested <- unique(gene_sets[[gene_set]])
  selected <- raw[gene %in% requested]
  selected[, gene_set := gene_set]
  selected[, program_order := match(gene, requested)]
  setorder(selected, program_order)
  gene_rows[[gene_set]] <- selected[, c(
    "gene_set", "program_order", "gene", effect_columns
  ), with = FALSE]

  for (block in blocks) {
    column <- paste0(block, "_log2_r_over_d")
    values <- selected[[column]]
    values <- values[is.finite(values)]
    summary_rows[[paste(gene_set, block, sep = "__")]] <- data.table(
      gene_set = gene_set,
      model = block,
      n_defined = length(requested),
      n_measurable = length(values),
      mean_log2_r_over_d = mean(values),
      median_log2_r_over_d = stats::median(values),
      q25_log2_r_over_d = as.numeric(stats::quantile(values, 0.25, names = FALSE)),
      q75_log2_r_over_d = as.numeric(stats::quantile(values, 0.75, names = FALSE)),
      positive_gene_n = sum(values > 0),
      positive_gene_fraction = mean(values > 0)
    )
  }
}

gene_effects <- rbindlist(gene_rows, use.names = TRUE)
gene_effects[, positive_all_four_models := rowSums(.SD > 0) == length(effect_columns),
             .SDcols = effect_columns]
model_summary <- rbindlist(summary_rows, use.names = TRUE)
all_four_summary <- gene_effects[, .(
  n_measurable = .N,
  positive_all_four_n = sum(positive_all_four_models),
  positive_all_four_fraction = mean(positive_all_four_models)
), by = gene_set]

write_tsv(gene_effects, "gene_effects.tsv")
write_tsv(model_summary, "model_summary.tsv")
write_tsv(all_four_summary, "all_four_model_summary.tsv")

run_info <- data.table(
  analysis_date = "2026-08-23",
  r_version = R.version.string,
  readxl_version = as.character(packageVersion("readxl")),
  data_table_version = as.character(packageVersion("data.table")),
  source_unit = "published_gene_level_model_summary",
  source_cohort = "57_pure_HGP_specimens_from_51_patients",
  workbook_md5 = unname(tools::md5sum(workbook)),
  program_md5 = unname(tools::md5sum(program_path)),
  core_md5 = unname(tools::md5sum(core_path)),
  spec_md5 = unname(tools::md5sum(spec_path))
)
write_tsv(run_info, "run_info.tsv")
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))

rec_uncorrected <- model_summary[
  gene_set == "frozen_rec_top50" & model == "uncorrected"
]
core_uncorrected <- model_summary[
  gene_set == "cross_modal_hgp_top15" & model == "uncorrected"
]
rec_all <- all_four_summary[gene_set == "frozen_rec_top50"]
summary_lines <- c(
  "# Fixed REC-program summary in the Latacz interface cohort",
  "",
  sprintf(
    "All %d/%d frozen REC genes were measurable.",
    rec_uncorrected$n_measurable, rec_uncorrected$n_defined
  ),
  "",
  "## Published uncorrected model",
  "",
  sprintf(
    paste0(
      "The mean gene effect was %+.3f log2(rHGP/dHGP), the median was %+.3f, ",
      "and %d/%d genes were positive."
    ),
    rec_uncorrected$mean_log2_r_over_d,
    rec_uncorrected$median_log2_r_over_d,
    rec_uncorrected$positive_gene_n,
    rec_uncorrected$n_measurable
  ),
  sprintf(
    paste0(
      "For the 15-gene cross-modal HGP core, the mean effect was %+.3f, the median was %+.3f, ",
      "and %d/%d genes were positive."
    ),
    core_uncorrected$mean_log2_r_over_d,
    core_uncorrected$median_log2_r_over_d,
    core_uncorrected$positive_gene_n,
    core_uncorrected$n_measurable
  ),
  "",
  "## Adjustment consistency",
  "",
  sprintf(
    "%d/%d REC genes remained positive in all four published models.",
    rec_all$positive_all_four_n, rec_all$n_measurable
  )
)
for (index in seq_len(nrow(model_summary[gene_set == "frozen_rec_top50"]))) {
  row <- model_summary[gene_set == "frozen_rec_top50"][index]
  summary_lines <- c(
    summary_lines,
    sprintf(
      "- %s: mean %+.3f, median %+.3f, %d/%d positive.",
      row$model, row$mean_log2_r_over_d, row$median_log2_r_over_d,
      row$positive_gene_n, row$n_measurable
    )
  )
}
summary_lines <- c(
  summary_lines,
  "",
  "## Interpretation boundary",
  "",
  paste(
    "These are summaries of correlated gene effects from the authors' models.",
    "They corroborate program direction at the tumour–liver interface but do not supply a patient-level program P value."
  )
)
writeLines(summary_lines, file.path(out_dir, "analysis_summary.md"), useBytes = TRUE)

output_manifest <- data.table(
  file = c(
    "gene_effects.tsv", "model_summary.tsv", "all_four_model_summary.tsv",
    "run_info.tsv", "sessionInfo.txt", "analysis_summary.md"
  ),
  role = c(
    "gene-level published effects", "fixed gene-set summaries by model",
    "four-model direction consistency", "source fingerprints and environment",
    "R environment", "human-readable result summary"
  )
)
fwrite(output_manifest, file.path(out_dir, "_analysis_outputs.md"), sep = "\t")
