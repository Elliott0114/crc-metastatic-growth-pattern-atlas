#!/usr/bin/env Rscript

# Analysis: External bulk and interface projection of frozen regulatory modules
# Date: 2026-09-01
# Random seed: 42
# Specification: metadata/external_regulatory_projection_spec_2026-09-01.md

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(readxl)
})

options(stringsAsFactors = FALSE)
random_seed <- 42L
bootstrap_iterations <- 10000L
set.seed(random_seed)

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
root <- if (length(file_arg)) {
  normalizePath(file.path(dirname(sub("^--file=", "", file_arg[[1]])), ".."), mustWork = TRUE)
} else {
  normalizePath(".", mustWork = TRUE)
}

out_dir <- file.path(
  root, "analysis_results", "deep_biology_upgrade_2026-08-31",
  "phase3_external_regulatory_projection"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

gene_set_path <- file.path(out_dir, "frozen_external_regulatory_gene_sets_long.tsv")
registry_path <- file.path(out_dir, "frozen_external_regulatory_gene_set_registry.tsv")
named_path <- file.path(out_dir, "frozen_external_regulatory_named_genes.tsv")
phase2_gene_set_path <- file.path(
  root, "analysis_results", "deep_biology_upgrade_2026-08-31",
  "phase2_mechanistic_specificity", "frozen_gene_sets_long.tsv"
)
bulk_path <- file.path(
  root, "data_sources", "GSE151165", "GSE151165_RNA_seq.raw_read_count.xlsx"
)
latacz_path <- file.path(
  root, "data_sources", "Latacz_2024_targeted_HGP_interface",
  "SupplTablesLataczforupload.xlsx"
)
spec_path <- file.path(
  root, "metadata", "external_regulatory_projection_spec_2026-09-01.md"
)
required_paths <- c(
  gene_set_path, registry_path, named_path, phase2_gene_set_path,
  bulk_path, latacz_path, spec_path
)
if (any(!file.exists(required_paths))) {
  stop("One or more required external-projection inputs are missing")
}

tol <- sqrt(.Machine$double.eps)

write_tsv <- function(x, filename) {
  fwrite(as.data.table(x), file.path(out_dir, filename), sep = "\t", na = "NA")
}

safe_gene_z <- function(x) {
  means <- rowMeans(x, na.rm = TRUE)
  sds <- apply(x, 1L, sd, na.rm = TRUE)
  keep <- is.finite(sds) & sds > 0
  z <- sweep(x[keep, , drop = FALSE], 1L, means[keep], "-")
  sweep(z, 1L, sds[keep], "/")
}

score_program <- function(gene_z, genes) {
  declared <- unique(toupper(as.character(genes)))
  measured <- intersect(declared, rownames(gene_z))
  list(
    score = if (length(measured)) {
      colMeans(gene_z[measured, , drop = FALSE], na.rm = TRUE)
    } else {
      NULL
    },
    declared = declared,
    measured = measured,
    missing = setdiff(declared, measured)
  )
}

hedges_g <- function(values, group) {
  r_values <- values[group == 1L]
  d_values <- values[group == 0L]
  n_r <- length(r_values)
  n_d <- length(d_values)
  pooled_variance <- (
    (n_r - 1) * var(r_values) + (n_d - 1) * var(d_values)
  ) / (n_r + n_d - 2)
  if (!is.finite(pooled_variance) || pooled_variance <= 0) return(NA_real_)
  correction <- 1 - 3 / (4 * (n_r + n_d) - 9)
  correction * (mean(r_values) - mean(d_values)) / sqrt(pooled_variance)
}

exact_two_group <- function(values, group) {
  values <- as.numeric(values)
  group <- as.integer(group)
  n_r <- sum(group == 1L)
  observed <- mean(values[group == 1L]) - mean(values[group == 0L])
  allocations <- combn(seq_along(values), n_r)
  r_sums <- colSums(matrix(values[allocations], nrow = n_r))
  null <- r_sums / n_r -
    (sum(values) - r_sums) / (length(values) - n_r)
  list(
    effect = observed,
    p_two_sided = mean(abs(null) >= abs(observed) - tol),
    p_one_sided_r_greater = mean(null >= observed - tol),
    n_allocations = length(null)
  )
}

bootstrap_group_difference <- function(values, group, iterations) {
  r_index <- which(group == 1L)
  d_index <- which(group == 0L)
  estimates <- replicate(iterations, {
    mean(values[sample(r_index, length(r_index), replace = TRUE)]) -
      mean(values[sample(d_index, length(d_index), replace = TRUE)])
  })
  as.numeric(quantile(estimates, c(0.025, 0.975), names = FALSE, type = 7))
}

group_coefficient <- function(outcome, group, covariate) {
  fit <- lm.fit(
    cbind(1, as.numeric(group), as.numeric(covariate)),
    as.numeric(outcome)
  )
  unname(fit$coefficients[[2L]])
}

exact_adjusted <- function(outcome, covariate, group) {
  outcome <- as.numeric(outcome)
  covariate <- as.numeric(covariate)
  group <- as.integer(group)
  observed <- group_coefficient(outcome, group, covariate)
  allocations <- combn(seq_along(group), sum(group == 1L))
  group_matrix <- matrix(0, nrow = length(group), ncol = ncol(allocations))
  group_matrix[cbind(
    as.vector(allocations),
    rep(seq_len(ncol(allocations)), each = nrow(allocations))
  )] <- 1
  nuisance <- cbind(1, covariate)
  residual_group <- group_matrix - nuisance %*%
    solve(crossprod(nuisance), crossprod(nuisance, group_matrix))
  denominators <- colSums(residual_group^2)
  null <- colSums(residual_group * outcome) / denominators
  null <- null[is.finite(null)]
  list(
    effect = observed,
    p_two_sided = mean(abs(null) >= abs(observed) - tol),
    n_allocations = length(null)
  )
}

bootstrap_adjusted <- function(outcome, covariate, group, iterations) {
  r_index <- which(group == 1L)
  d_index <- which(group == 0L)
  estimates <- replicate(iterations, {
    index <- c(
      sample(r_index, length(r_index), replace = TRUE),
      sample(d_index, length(d_index), replace = TRUE)
    )
    group_coefficient(outcome[index], group[index], covariate[index])
  })
  estimates <- estimates[is.finite(estimates)]
  if (length(estimates) < iterations * 0.95) return(c(NA_real_, NA_real_))
  as.numeric(quantile(estimates, c(0.025, 0.975), names = FALSE, type = 7))
}

gene_sets_long <- fread(gene_set_path)
registry <- fread(registry_path)
named_genes <- fread(named_path)
phase2_gene_sets_long <- fread(phase2_gene_set_path)
gene_sets_long[, gene := toupper(gene)]
named_genes[, gene := toupper(gene)]
phase2_gene_sets_long[, gene := toupper(gene)]
gene_sets <- split(gene_sets_long$gene, gene_sets_long$set_id)
gene_set_order <- registry$set_id
if (!identical(names(gene_sets)[match(gene_set_order, names(gene_sets))], gene_set_order)) {
  stop("Frozen gene-set registry and membership table do not align")
}

composition_sets <- split(
  phase2_gene_sets_long[set_id %chin% c("PAN_EPITHELIAL_7", "HEPATOCYTE_CONTEXT_6"), gene],
  phase2_gene_sets_long[set_id %chin% c("PAN_EPITHELIAL_7", "HEPATOCYTE_CONTEXT_6"), set_id]
)
if (!all(c("PAN_EPITHELIAL_7", "HEPATOCYTE_CONTEXT_6") %in% names(composition_sets))) {
  stop("Frozen composition-control panels are missing")
}

# GSE151165 tumour-level TMM log2-CPM.
raw_bulk <- as.data.table(read_excel(bulk_path, .name_repair = "unique_quiet"))
bulk_samples <- grep("^KR-", names(raw_bulk), value = TRUE)
if (length(bulk_samples) != 30L) stop("Expected 30 GSE151165 sample columns")
setnames(raw_bulk, 1L, "source_gene")
raw_bulk[, gene := toupper(trimws(as.character(source_gene)))]
raw_bulk <- raw_bulk[!is.na(gene) & gene != ""]
bulk_counts_dt <- raw_bulk[
  , lapply(.SD, function(x) sum(as.numeric(x), na.rm = TRUE)),
  by = gene,
  .SDcols = bulk_samples
]
bulk_counts <- as.matrix(bulk_counts_dt[, ..bulk_samples])
storage.mode(bulk_counts) <- "double"
rownames(bulk_counts) <- bulk_counts_dt$gene

bulk_metadata_all <- data.table(
  sample_id = bulk_samples,
  source_class = c(rep("D_N", 9L), rep("D_T", 9L), rep("R_N", 6L), rep("R_T", 6L))
)
bulk_metadata_all[, `:=`(
  hgp = fifelse(startsWith(source_class, "R"), "rHGP", "dHGP"),
  tissue = fifelse(endsWith(source_class, "_T"), "tumour", "adjacent_liver")
)]
bulk_metadata <- bulk_metadata_all[tissue == "tumour"]
bulk_counts_tumour <- bulk_counts[, bulk_metadata$sample_id, drop = FALSE]
bulk_dge <- DGEList(bulk_counts_tumour, group = bulk_metadata$hgp)
bulk_keep <- filterByExpr(bulk_dge, group = bulk_metadata$hgp)
bulk_dge <- calcNormFactors(bulk_dge[bulk_keep, , keep.lib.sizes = FALSE])
effective_library <- bulk_dge$samples$lib.size * bulk_dge$samples$norm.factors
names(effective_library) <- rownames(bulk_dge$samples)
bulk_log_cpm <- cpm(
  bulk_counts_tumour,
  lib.size = effective_library[colnames(bulk_counts_tumour)],
  log = TRUE,
  prior.count = 1
)
bulk_gene_z <- safe_gene_z(bulk_log_cpm)
bulk_group <- as.integer(bulk_metadata$hgp == "rHGP")

score_rows <- list()
coverage_rows <- list()
contrast_rows <- list()
loo_rows <- list()
module_scores <- list()
for (set_id in gene_set_order) {
  scored <- score_program(bulk_gene_z, gene_sets[[set_id]])
  if (is.null(scored$score)) stop(paste("No measurable genes for", set_id))
  values <- as.numeric(scored$score[bulk_metadata$sample_id])
  module_scores[[set_id]] <- values
  coverage_rows[[set_id]] <- data.table(
    set_id = set_id,
    n_defined = length(scored$declared),
    n_measured_nonconstant = length(scored$measured),
    coverage_fraction = length(scored$measured) / length(scored$declared),
    measured_genes = paste(scored$measured, collapse = ";"),
    missing_or_zero_variance_genes = paste(scored$missing, collapse = ";")
  )
  score_rows[[set_id]] <- data.table(
    sample_id = bulk_metadata$sample_id,
    hgp = bulk_metadata$hgp,
    set_id = set_id,
    score = values
  )
  exact <- exact_two_group(values, bulk_group)
  interval <- bootstrap_group_difference(values, bulk_group, bootstrap_iterations)
  contrast_rows[[set_id]] <- data.table(
    set_id = set_id,
    n_patients = length(values),
    n_rhgp = sum(bulk_group == 1L),
    n_dhgp = sum(bulk_group == 0L),
    rhgp_mean = mean(values[bulk_group == 1L]),
    dhgp_mean = mean(values[bulk_group == 0L]),
    rhgp_minus_dhgp = exact$effect,
    bootstrap_ci_lower = interval[[1L]],
    bootstrap_ci_upper = interval[[2L]],
    hedges_g = hedges_g(values, bulk_group),
    exact_p_one_sided_r_greater = exact$p_one_sided_r_greater,
    exact_p_two_sided = exact$p_two_sided,
    n_exact_allocations = exact$n_allocations
  )
  loo_rows[[set_id]] <- rbindlist(lapply(seq_along(values), function(index) {
    retained <- seq_along(values) != index
    data.table(
      set_id = set_id,
      omitted_sample = bulk_metadata$sample_id[[index]],
      omitted_hgp = bulk_metadata$hgp[[index]],
      rhgp_minus_dhgp = mean(values[retained & bulk_group == 1L]) -
        mean(values[retained & bulk_group == 0L])
    )
  }))
}

bulk_scores <- rbindlist(score_rows, use.names = TRUE)
bulk_coverage <- rbindlist(coverage_rows, use.names = TRUE)
bulk_contrasts <- rbindlist(contrast_rows, use.names = TRUE)
bulk_contrasts[, fdr_bh_all_modules := p.adjust(exact_p_two_sided, method = "BH")]
bulk_loo <- rbindlist(loo_rows, use.names = TRUE)
bulk_loo[, positive_direction := rhgp_minus_dhgp > 0]

composition_scores <- lapply(composition_sets, function(genes) {
  scored <- score_program(bulk_gene_z, genes)
  if (is.null(scored$score)) stop("Composition panel could not be scored")
  as.numeric(scored$score[bulk_metadata$sample_id])
})
adjusted_rows <- list()
row_index <- 0L
for (covariate_id in names(composition_scores)) {
  covariate <- composition_scores[[covariate_id]]
  for (set_id in gene_set_order) {
    row_index <- row_index + 1L
    outcome <- module_scores[[set_id]]
    exact <- exact_adjusted(outcome, covariate, bulk_group)
    interval <- bootstrap_adjusted(
      outcome, covariate, bulk_group, bootstrap_iterations
    )
    adjusted_rows[[row_index]] <- data.table(
      set_id = set_id,
      adjustment = covariate_id,
      adjusted_hgp_coefficient = exact$effect,
      bootstrap_ci_lower = interval[[1L]],
      bootstrap_ci_upper = interval[[2L]],
      exact_p_two_sided = exact$p_two_sided,
      n_exact_allocations = exact$n_allocations,
      hgp_covariate_correlation = cor(
        bulk_group, covariate, method = "spearman"
      ),
      outcome_covariate_correlation = cor(
        outcome, covariate, method = "spearman"
      )
    )
  }
}
bulk_adjusted <- rbindlist(adjusted_rows, use.names = TRUE)
bulk_adjusted[, fdr_bh_within_adjustment := p.adjust(exact_p_two_sided, "BH"), by = adjustment]

named_contrast_rows <- list()
named_value_rows <- list()
for (index in seq_len(nrow(named_genes))) {
  gene <- named_genes$gene[[index]]
  if (!gene %in% rownames(bulk_log_cpm)) next
  values <- as.numeric(bulk_log_cpm[gene, bulk_metadata$sample_id])
  exact <- exact_two_group(values, bulk_group)
  interval <- bootstrap_group_difference(values, bulk_group, bootstrap_iterations)
  raw_values <- as.numeric(bulk_counts_tumour[gene, bulk_metadata$sample_id])
  named_contrast_rows[[gene]] <- data.table(
    gene = gene,
    role = named_genes$role[[index]],
    rhgp_mean_log2_cpm = mean(values[bulk_group == 1L]),
    dhgp_mean_log2_cpm = mean(values[bulk_group == 0L]),
    rhgp_minus_dhgp_log2_cpm = exact$effect,
    bootstrap_ci_lower = interval[[1L]],
    bootstrap_ci_upper = interval[[2L]],
    hedges_g = hedges_g(values, bulk_group),
    exact_p_two_sided = exact$p_two_sided,
    n_exact_allocations = exact$n_allocations,
    detected_rhgp_n = sum(raw_values[bulk_group == 1L] > 0),
    detected_dhgp_n = sum(raw_values[bulk_group == 0L] > 0)
  )
  named_value_rows[[gene]] <- data.table(
    sample_id = bulk_metadata$sample_id,
    hgp = bulk_metadata$hgp,
    gene = gene,
    raw_count = raw_values,
    log2_cpm = values
  )
}
bulk_named_contrasts <- rbindlist(named_contrast_rows, use.names = TRUE)
bulk_named_contrasts[, fdr_bh_named_genes := p.adjust(exact_p_two_sided, "BH")]
bulk_named_values <- rbindlist(named_value_rows, use.names = TRUE)


write_tsv(bulk_contrasts, "gse151165_regulatory_module_contrasts.tsv")
