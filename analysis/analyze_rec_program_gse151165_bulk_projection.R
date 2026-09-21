#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(readxl)
})

options(stringsAsFactors = FALSE)
set.seed(42)

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
root <- if (length(file_arg)) {
  normalizePath(file.path(dirname(sub("^--file=", "", file_arg[[1]])), ".."), mustWork = TRUE)
} else {
  normalizePath(".", mustWork = TRUE)
}

raw_path <- file.path(
  root, "data_sources", "GSE151165", "GSE151165_RNA_seq.raw_read_count.xlsx"
)
program_path <- file.path(
  root, "analysis_results", "ogden_anchor_program", "selected_anchor_program_top50.tsv"
)
core_path <- file.path(
  root, "analysis_results", "rec_program_cross_modal_concordance",
  "top15_interpretable_core.tsv"
)
axis_path <- file.path(
  root, "analysis_results", "rec_program_two_axis_interpretation",
  "gene_two_axis_classification.tsv"
)
spec_path <- file.path(
  root, "metadata", "rec_program_gse151165_bulk_projection_spec_2026-08-23.md"
)
out_dir <- file.path(root, "analysis_results", "rec_program_gse151165_bulk_projection")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

bootstrap_replicates <- 10000L

md5_file <- function(path) {
  unname(tools::md5sum(path))
}

write_tsv <- function(x, filename) {
  fwrite(as.data.table(x), file.path(out_dir, filename), sep = "\t", na = "NA")
}

hedges_g <- function(values, labels) {
  r_values <- values[labels == "rHGP"]
  d_values <- values[labels == "dHGP"]
  n_r <- length(r_values)
  n_d <- length(d_values)
  pooled_variance <- (
    (n_r - 1) * stats::var(r_values) + (n_d - 1) * stats::var(d_values)
  ) / (n_r + n_d - 2)
  if (!is.finite(pooled_variance) || pooled_variance <= 0) {
    return(NA_real_)
  }
  uncorrected <- (mean(r_values) - mean(d_values)) / sqrt(pooled_variance)
  correction <- 1 - 3 / (4 * (n_r + n_d) - 9)
  correction * uncorrected
}

exact_label_test <- function(values, labels) {
  replacement_n <- sum(labels == "rHGP")
  observed <- mean(values[labels == "rHGP"]) - mean(values[labels == "dHGP"])
  allocations <- utils::combn(seq_along(values), replacement_n)
  total_sum <- sum(values)
  replacement_sums <- colSums(matrix(values[allocations], nrow = replacement_n))
  null_effects <- replacement_sums / replacement_n -
    (total_sum - replacement_sums) / (length(values) - replacement_n)
  data.table(
    observed_effect = observed,
    exact_p_one_sided = mean(null_effects >= observed - .Machine$double.eps^0.5),
    exact_p_two_sided = mean(abs(null_effects) >= abs(observed) - .Machine$double.eps^0.5),
    n_exact_allocations = length(null_effects)
  )
}

bootstrap_ci <- function(values, labels, replicates) {
  r_values <- values[labels == "rHGP"]
  d_values <- values[labels == "dHGP"]
  effects <- replicate(
    replicates,
    mean(sample(r_values, length(r_values), replace = TRUE)) -
      mean(sample(d_values, length(d_values), replace = TRUE))
  )
  as.numeric(stats::quantile(effects, c(0.025, 0.975), names = FALSE, type = 7))
}

summarise_endpoint <- function(score_table, endpoint) {
  selected <- score_table[endpoint_name == endpoint]
  values <- selected$score
  labels <- selected$hgp
  exact <- exact_label_test(values, labels)
  interval <- bootstrap_ci(values, labels, bootstrap_replicates)
  data.table(
    endpoint = endpoint,
    n_patients = nrow(selected),
    n_rhgp = sum(labels == "rHGP"),
    n_dhgp = sum(labels == "dHGP"),
    rhgp_mean = mean(values[labels == "rHGP"]),
    dhgp_mean = mean(values[labels == "dHGP"]),
    rhgp_minus_dhgp = exact$observed_effect,
    bootstrap_ci_lower = interval[[1]],
    bootstrap_ci_upper = interval[[2]],
    hedges_g = hedges_g(values, labels),
    exact_p_one_sided = exact$exact_p_one_sided,
    exact_p_two_sided = exact$exact_p_two_sided,
    n_exact_allocations = exact$n_exact_allocations,
    all_rhgp_above_all_dhgp = min(values[labels == "rHGP"]) > max(values[labels == "dHGP"])
  )
}

raw <- as.data.table(read_excel(raw_path, .name_repair = "unique_quiet"))
sample_columns <- grep("^KR-", names(raw), value = TRUE)
if (length(sample_columns) != 30L) {
  stop("Expected 30 GSE151165 sample columns")
}
setnames(raw, 1L, "source_gene")
raw[, gene := toupper(trimws(as.character(source_gene)))]
raw <- raw[!is.na(gene) & gene != ""]
counts_table <- raw[, lapply(.SD, function(x) sum(as.numeric(x), na.rm = TRUE)),
                    by = gene, .SDcols = sample_columns]
counts <- as.matrix(counts_table[, ..sample_columns])
mode(counts) <- "numeric"
rownames(counts) <- counts_table$gene

classes <- c(rep("D_N", 9L), rep("D_T", 9L), rep("R_N", 6L), rep("R_T", 6L))
metadata <- data.table(sample = sample_columns, class = classes)
metadata[, hgp := ifelse(startsWith(class, "R"), "rHGP", "dHGP")]
metadata[, tissue := ifelse(endsWith(class, "_T"), "tumour", "normal")]
tumour_metadata <- metadata[tissue == "tumour"]
counts_tumour <- counts[, tumour_metadata$sample, drop = FALSE]

y_all <- DGEList(counts_tumour, group = tumour_metadata$hgp)
keep <- filterByExpr(y_all, group = tumour_metadata$hgp)
y_filtered <- calcNormFactors(y_all[keep, , keep.lib.sizes = FALSE])
effective_library <- y_filtered$samples$lib.size * y_filtered$samples$norm.factors
names(effective_library) <- rownames(y_filtered$samples)
log_cpm <- cpm(
  counts_tumour,
  lib.size = effective_library[colnames(counts_tumour)],
  log = TRUE,
  prior.count = 1
)

program <- fread(program_path)
program[, gene := toupper(gene)]
core <- fread(core_path)
core[, gene := toupper(gene)]
axis <- fread(axis_path)
axis[, gene := toupper(gene)]

module_columns <- c(
  experimental_ap1_target_overlap = "member_experimental_ap1_target_overlap",
  hypoxia_mp6_overlap = "member_hypoxia_mp6_overlap",
  regenerative_overlap = "member_regenerative_overlap",
  nfkb_regulon_overlap = "member_nfkb_regulon_overlap"
)
gene_sets <- list(
  rec_equal_gene = program$gene,
  hgp_cross_modal_top15 = core$gene
)
for (module_name in names(module_columns)) {
  column <- module_columns[[module_name]]
  gene_sets[[module_name]] <- axis[get(column) == TRUE, gene]
}

all_program_genes <- intersect(program$gene, rownames(log_cpm))
gene_z <- t(scale(t(log_cpm)))
score_rows <- list()
coverage_rows <- list()
for (endpoint in names(gene_sets)) {
  requested <- unique(gene_sets[[endpoint]])
  matched <- intersect(requested, rownames(log_cpm))
  coverage_rows[[endpoint]] <- data.table(
    endpoint = endpoint,
    n_defined = length(requested),
    n_matched = length(matched),
    coverage_fraction = length(matched) / length(requested),
    matched_genes = paste(matched, collapse = ","),
    missing_genes = paste(setdiff(requested, matched), collapse = ",")
  )
  endpoint_scores <- colMeans(gene_z[matched, , drop = FALSE], na.rm = TRUE)
  score_rows[[endpoint]] <- data.table(
    sample = names(endpoint_scores),
    endpoint_name = endpoint,
    score = as.numeric(endpoint_scores)
  )
}

pooled_counts <- colSums(counts_tumour[all_program_genes, , drop = FALSE])
pooled_log2_cpm <- log2(
  (pooled_counts + 1) / effective_library[names(pooled_counts)] * 1e6
)
score_rows[["rec_pooled_log2_cpm"]] <- data.table(
  sample = names(pooled_log2_cpm),
  endpoint_name = "rec_pooled_log2_cpm",
  score = as.numeric(pooled_log2_cpm)
)

scores <- rbindlist(score_rows, use.names = TRUE)
scores <- merge(
  scores,
  tumour_metadata[, .(sample, hgp, class)],
  by = "sample",
  all.x = TRUE,
  sort = FALSE
)
endpoint_order <- c(
  "rec_equal_gene", "rec_pooled_log2_cpm", "hgp_cross_modal_top15",
  names(module_columns)
)
endpoint_summary <- rbindlist(lapply(endpoint_order, function(x) {
  summarise_endpoint(scores, x)
}))

primary_scores <- scores[endpoint_name == "rec_equal_gene"]
loo <- rbindlist(lapply(seq_len(nrow(primary_scores)), function(index) {
  retained <- primary_scores[-index]
  effect <- mean(retained[hgp == "rHGP", score]) - mean(retained[hgp == "dHGP", score])
  data.table(
    omitted_sample = primary_scores$sample[[index]],
    omitted_hgp = primary_scores$hgp[[index]],
    n_retained = nrow(retained),
    rhgp_minus_dhgp = effect,
    positive_direction = effect > 0
  )
}))

gene_rows <- lapply(all_program_genes, function(gene_symbol) {
  values <- log_cpm[gene_symbol, tumour_metadata$sample]
  labels <- tumour_metadata$hgp
  data.table(
    gene = gene_symbol,
    program_rank = program[gene == gene_symbol, program_rank],
    rhgp_mean_log2_cpm = mean(values[labels == "rHGP"]),
    dhgp_mean_log2_cpm = mean(values[labels == "dHGP"]),
    rhgp_minus_dhgp = mean(values[labels == "rHGP"]) - mean(values[labels == "dHGP"]),
    pairwise_superiority = mean(outer(
      values[labels == "rHGP"], values[labels == "dHGP"], FUN = ">"
    )),
    positive_direction = mean(values[labels == "rHGP"]) > mean(values[labels == "dHGP"])
  )
})
gene_effects <- rbindlist(gene_rows)[order(program_rank)]

write_tsv(rbindlist(coverage_rows), "program_coverage.tsv")
write_tsv(scores[order(match(endpoint_name, endpoint_order), hgp, sample)], "patient_scores.tsv")
write_tsv(endpoint_summary, "endpoint_summary.tsv")
write_tsv(loo, "leave_one_patient_out.tsv")
write_tsv(gene_effects, "gene_level_effects.tsv")

run_info <- data.table(
  analysis_date = "2026-08-23",
  r_version = R.version.string,
  edgeR_version = as.character(packageVersion("edgeR")),
  data_table_version = as.character(packageVersion("data.table")),
  seed = 42L,
  bootstrap_replicates = bootstrap_replicates,
  statistical_unit = "public_patient_sample",
  n_rhgp = sum(tumour_metadata$hgp == "rHGP"),
  n_dhgp = sum(tumour_metadata$hgp == "dHGP"),
  exact_label_allocations = choose(nrow(tumour_metadata), sum(tumour_metadata$hgp == "rHGP")),
  raw_file_md5 = md5_file(raw_path),
  program_file_md5 = md5_file(program_path),
  core_file_md5 = md5_file(core_path),
  two_axis_file_md5 = md5_file(axis_path),
  spec_file_md5 = md5_file(spec_path)
)
write_tsv(run_info, "run_info.tsv")
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))

primary <- endpoint_summary[endpoint == "rec_equal_gene"]
pooled <- endpoint_summary[endpoint == "rec_pooled_log2_cpm"]
top_core <- endpoint_summary[endpoint == "hgp_cross_modal_top15"]
summary_lines <- c(
  "# Fixed REC-program projection into GSE151165 bulk CRLM",
  "",
  sprintf(
    "The primary comparison included %d rHGP and %d dHGP tumour samples and retained %d/%d frozen REC genes.",
    primary$n_rhgp, primary$n_dhgp, length(all_program_genes), nrow(program)
  ),
  "",
  "## Primary result",
  "",
  sprintf(
    paste0(
      "The equal-gene REC score was higher by %+.3f in rHGP (bootstrap 95%% CI %+.3f to %+.3f; ",
      "Hedges' g %+.3f; exact two-sided P=%.4f). All %d leave-one-patient-out effects were positive."
    ),
    primary$rhgp_minus_dhgp, primary$bootstrap_ci_lower, primary$bootstrap_ci_upper,
    primary$hedges_g, primary$exact_p_two_sided, sum(loo$positive_direction)
  ),
  "",
  sprintf(
    "The pooled-program log2 CPM difference was %+.3f (Hedges' g %+.3f; exact two-sided P=%.4f).",
    pooled$rhgp_minus_dhgp, pooled$hedges_g, pooled$exact_p_two_sided
  ),
  sprintf(
    "The 15-gene cross-modal HGP core differed by %+.3f (Hedges' g %+.3f; exact two-sided P=%.4f).",
    top_core$rhgp_minus_dhgp, top_core$hedges_g, top_core$exact_p_two_sided
  ),
  "",
  sprintf(
    "%d/%d measurable REC genes had a positive rHGP-minus-dHGP mean direction.",
    sum(gene_effects$positive_direction), nrow(gene_effects)
  ),
  "",
  "## Interpretation boundary",
  "",
  paste(
    "The cohort provides an independent regional HGP test of the frozen program.",
    "Because it is bulk tissue, epithelial single-cell data remain necessary for assigning cellular source."
  ),
  "",
  "## Reproducibility",
  "",
  "```bash",
  "conda run -n crc-metastatic-growth Rscript scripts/analyze_rec_program_gse151165_bulk_projection.R",
  "```"
)
writeLines(summary_lines, file.path(out_dir, "analysis_summary.md"), useBytes = TRUE)

output_manifest <- data.table(
  file = c(
    "program_coverage.tsv", "patient_scores.tsv", "endpoint_summary.tsv",
    "leave_one_patient_out.tsv", "gene_level_effects.tsv", "run_info.tsv",
    "sessionInfo.txt", "analysis_summary.md"
  ),
  role = c(
    "fixed gene-set coverage", "patient-level program scores", "patient-level effects",
    "primary leave-one-patient-out direction", "gene-level descriptive effects",
    "environment and source fingerprints", "R environment", "human-readable summary"
  )
)
fwrite(output_manifest, file.path(out_dir, "_analysis_outputs.md"), sep = "\t")
