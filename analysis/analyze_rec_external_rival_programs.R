#!/usr/bin/env Rscript

# Analysis: External one-rival-at-a-time controls for the frozen REC programme
# Date: 2026-09-01
# Random seed: 42
# R and package versions are written to sessionInfo output.

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(readxl)
})

random_seed <- 42L
set.seed(random_seed)

root <- normalizePath(getwd(), mustWork = TRUE)
out_dir <- file.path(
  root,
  "analysis_results",
  "deep_biology_upgrade_2026-08-31",
  "phase2_mechanistic_specificity"
)
registry_path <- file.path(out_dir, "frozen_gene_set_registry.tsv")
gene_set_path <- file.path(out_dir, "frozen_gene_sets_long.tsv")
bulk_path <- file.path(
  root, "data_sources", "GSE151165", "GSE151165_RNA_seq.raw_read_count.xlsx"
)
spatial_counts_path <- file.path(
  root, "analysis_results", "rec_program_reader_facing_specificity",
  "gse294385_patient_region_all_gene_counts.tsv.gz"
)
spatial_metadata_path <- file.path(
  root, "analysis_results", "rec_program_reader_facing_specificity",
  "gse294385_patient_region_metadata.tsv"
)
spatial_manifest_path <- file.path(
  root, "metadata", "gse294385_liver_sample_manifest.tsv"
)
required_paths <- c(
  registry_path, gene_set_path, bulk_path, spatial_counts_path,
  spatial_metadata_path, spatial_manifest_path
)
if (any(!file.exists(required_paths))) {
  stop("One or more external rival-program inputs are missing")
}

bootstrap_iterations <- 5000L
tol <- sqrt(.Machine$double.eps)

stable_model_seed <- function(endpoint, rival_id, analysis) {
  key <- paste(endpoint, rival_id, analysis, sep = "::")
  code_points <- utf8ToInt(key)
  offset <- sum(as.double(code_points) * (seq_along(code_points) + 97))
  as.integer((random_seed + offset) %% 2147483646 + 1)
}

safe_z_rows <- function(x, population_sd = FALSE) {
  means <- rowMeans(x, na.rm = TRUE)
  if (population_sd) {
    centered <- sweep(x, 1L, means, "-")
    sds <- sqrt(rowMeans(centered^2, na.rm = TRUE))
  } else {
    sds <- apply(x, 1L, sd, na.rm = TRUE)
  }
  keep <- is.finite(sds) & sds > 0
  z <- sweep(x[keep, , drop = FALSE], 1L, means[keep], "-")
  sweep(z, 1L, sds[keep], "/")
}

score_program <- function(gene_z, genes) {
  declared <- unique(toupper(as.character(genes)))
  matched <- intersect(declared, rownames(gene_z))
  if (length(matched) < 3L) {
    return(list(score = NULL, n_declared = length(declared), n_measured = length(matched)))
  }
  list(
    score = colMeans(gene_z[matched, , drop = FALSE], na.rm = TRUE),
    n_declared = length(declared),
    n_measured = length(matched)
  )
}

exact_two_group <- function(values, group) {
  values <- as.numeric(values)
  group <- as.integer(group)
  n_positive <- sum(group == 1L)
  observed <- mean(values[group == 1L]) - mean(values[group == 0L])
  allocations <- combn(seq_along(values), n_positive)
  positive_sums <- colSums(matrix(values[allocations], nrow = n_positive))
  null <- positive_sums / n_positive -
    (sum(values) - positive_sums) / (length(values) - n_positive)
  list(
    effect = observed,
    exact_p_two_sided = mean(abs(null) >= abs(observed) - tol),
    n_assignments = length(null)
  )
}

exact_paired <- function(differences) {
  differences <- as.numeric(differences)
  signs <- as.matrix(expand.grid(rep(list(c(-1, 1)), length(differences))))
  null <- rowMeans(sweep(signs, 2L, differences, "*"))
  observed <- mean(differences)
  list(
    effect = observed,
    exact_p_two_sided = mean(abs(null) >= abs(observed) - tol),
    n_assignments = length(null)
  )
}

group_coefficient <- function(outcome, group, covariate) {
  fit <- lm.fit(cbind(1, as.numeric(group), as.numeric(covariate)), as.numeric(outcome))
  unname(fit$coefficients[[2L]])
}

intercept_coefficient <- function(outcome_difference, covariate_difference) {
  fit <- lm.fit(
    cbind(1, as.numeric(covariate_difference)),
    as.numeric(outcome_difference)
  )
  unname(fit$coefficients[[1L]])
}

exact_adjusted_unpaired <- function(outcome, covariate, group) {
  outcome <- as.numeric(outcome)
  covariate <- as.numeric(covariate)
  group <- as.integer(group)
  observed <- group_coefficient(outcome, group, covariate)
  allocations <- combn(seq_along(group), sum(group == 1L))
  group_matrix <- matrix(0, nrow = length(group), ncol = ncol(allocations))
  group_matrix[cbind(as.vector(allocations), rep(seq_len(ncol(allocations)), each = nrow(allocations)))] <- 1
  nuisance <- cbind(1, covariate)
  residual_group <- group_matrix - nuisance %*%
    solve(crossprod(nuisance), crossprod(nuisance, group_matrix))
  denominators <- colSums(residual_group^2)
  null <- colSums(residual_group * outcome) / denominators
  null <- null[is.finite(null)]
  list(
    effect = observed,
    exact_p_two_sided = mean(abs(null) >= abs(observed) - tol),
    n_assignments = length(null)
  )
}

exact_adjusted_paired <- function(outcome_difference, covariate_difference) {
  outcome_difference <- as.numeric(outcome_difference)
  covariate_difference <- as.numeric(covariate_difference)
  observed <- intercept_coefficient(outcome_difference, covariate_difference)
  signs <- as.matrix(expand.grid(rep(list(c(-1, 1)), length(outcome_difference))))
  null <- apply(signs, 1L, function(sign_vector) {
    intercept_coefficient(
      outcome_difference * sign_vector,
      covariate_difference * sign_vector
    )
  })
  null <- null[is.finite(null)]
  list(
    effect = observed,
    exact_p_two_sided = mean(abs(null) >= abs(observed) - tol),
    n_assignments = length(null)
  )
}

bootstrap_unpaired_ci <- function(outcome, covariate, group, iterations) {
  positive <- which(group == 1L)
  negative <- which(group == 0L)
  estimates <- replicate(iterations, {
    index <- c(
      sample(positive, length(positive), replace = TRUE),
      sample(negative, length(negative), replace = TRUE)
    )
    group_coefficient(outcome[index], group[index], covariate[index])
  })
  estimates <- estimates[is.finite(estimates)]
  if (length(estimates) < iterations * 0.95) return(c(NA_real_, NA_real_))
  quantile(estimates, c(0.025, 0.975), names = FALSE, type = 6)
}

bootstrap_paired_ci <- function(outcome_difference, covariate_difference, iterations) {
  n <- length(outcome_difference)
  estimates <- replicate(iterations, {
    index <- sample.int(n, n, replace = TRUE)
    intercept_coefficient(outcome_difference[index], covariate_difference[index])
  })
  estimates <- estimates[is.finite(estimates)]
  if (length(estimates) < iterations * 0.95) return(c(NA_real_, NA_real_))
  quantile(estimates, c(0.025, 0.975), names = FALSE, type = 6)
}

bootstrap_mean_ci <- function(values, iterations) {
  n <- length(values)
  estimates <- replicate(iterations, mean(sample(values, n, replace = TRUE)))
  quantile(estimates, c(0.025, 0.975), names = FALSE, type = 6)
}

registry <- fread(registry_path)
gene_sets_long <- fread(gene_set_path)
gene_sets <- split(toupper(gene_sets_long$gene), gene_sets_long$set_id)
if (!"FROZEN_REC_TOP50" %in% names(gene_sets)) {
  stop("Frozen REC programme is missing from the registry")
}
rec_genes <- unique(gene_sets[["FROZEN_REC_TOP50"]])
rival_registry <- registry[
  set_id != "FROZEN_REC_TOP50" & role != "overlap_audit_only"
]

# GSE151165 tumour-level TMM log-CPM.
raw_bulk <- as.data.table(read_excel(bulk_path, .name_repair = "unique_quiet"))
bulk_samples <- grep("^KR-", names(raw_bulk), value = TRUE)
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
bulk_effective_library <- bulk_dge$samples$lib.size * bulk_dge$samples$norm.factors
names(bulk_effective_library) <- rownames(bulk_dge$samples)
bulk_log_cpm <- cpm(
  bulk_counts_tumour,
  lib.size = bulk_effective_library[colnames(bulk_counts_tumour)],
  log = TRUE,
  prior.count = 1
)
bulk_gene_z <- safe_z_rows(bulk_log_cpm, population_sd = FALSE)

# GSE294385 patient-region pseudobulks, matching the audited existing projection.
spatial_long <- fread(cmd = paste("gzip -dc", shQuote(spatial_counts_path)))
spatial_metadata <- fread(spatial_metadata_path)
spatial_long[, sample_id := paste(patient, region, sep = "__")]
spatial_metadata[, sample_id := paste(patient, region, sep = "__")]
spatial_wide <- dcast(spatial_long, gene ~ sample_id, value.var = "raw_count", fill = 0)
spatial_genes <- toupper(spatial_wide$gene)
spatial_counts <- as.matrix(spatial_wide[, -"gene"])
storage.mode(spatial_counts) <- "double"
rownames(spatial_counts) <- spatial_genes
spatial_metadata <- spatial_metadata[match(colnames(spatial_counts), sample_id)]
if (anyNA(spatial_metadata$sample_id)) stop("Spatial metadata ordering failed")

spatial_manifest <- fread(spatial_manifest_path)
spatial_manifest <- spatial_manifest[
  tolower(as.character(selected_for_paired_extension)) == "true"
]
spatial_feature_sets <- lapply(spatial_manifest$sample, function(sample_id) {
  feature_path <- file.path(
    root, "data_sources", "Liu_2026_GSE294385", "extracted",
    sample_id, "filtered_feature_bc_matrix", "features.tsv.gz"
  )
  unique(toupper(fread(
    cmd = paste("gzip -dc", shQuote(feature_path)),
    header = FALSE,
    select = 2L
  )[[1L]]))
})
common_spatial_features <- Reduce(intersect, spatial_feature_sets)
spatial_counts <- spatial_counts[
  intersect(rownames(spatial_counts), common_spatial_features),,
  drop = FALSE
]
spatial_log_cpm <- log2(
  sweep(
    spatial_counts + 0.5,
    2L,
    spatial_metadata$library_size + 1,
    "/"
  ) * 1e6
)
spatial_gene_z <- safe_z_rows(spatial_log_cpm, population_sd = TRUE)

bulk_rec <- score_program(bulk_gene_z, rec_genes)
spatial_rec <- score_program(spatial_gene_z, rec_genes)
if (is.null(bulk_rec$score) || is.null(spatial_rec$score)) {
  stop("Frozen REC programme could not be scored in both external datasets")
}

bulk_group <- as.integer(bulk_metadata$hgp == "rHGP")
bulk_rec_score <- as.numeric(bulk_rec$score[bulk_metadata$sample_id])
bulk_baseline <- exact_two_group(bulk_rec_score, bulk_group)
bulk_baseline_ci <- c(
  quantile(
    replicate(bootstrap_iterations, {
      positive <- sample(bulk_rec_score[bulk_group == 1L], sum(bulk_group == 1L), TRUE)
      negative <- sample(bulk_rec_score[bulk_group == 0L], sum(bulk_group == 0L), TRUE)
      mean(positive) - mean(negative)
    }),
    c(0.025, 0.975), names = FALSE, type = 6
  )
)

spatial_score_dt <- copy(spatial_metadata)
spatial_score_dt[, rec_score := as.numeric(spatial_rec$score[sample_id])]
spatial_rec_wide <- dcast(spatial_score_dt, patient ~ region, value.var = "rec_score")
spatial_rec_wide[, rec_difference := macro_tumour - micro_tumour]
spatial_baseline <- exact_paired(spatial_rec_wide$rec_difference)
spatial_baseline_ci <- bootstrap_mean_ci(
  spatial_rec_wide$rec_difference,
  bootstrap_iterations
)

baseline <- rbind(
  data.table(
    dataset = "GSE151165 HGP-labelled bulk tumour",
    endpoint = "rHGP_minus_dHGP",
    n_units = nrow(bulk_metadata),
    n_rec_genes_defined = length(rec_genes),
    n_rec_genes_measured = bulk_rec$n_measured,
    rec_effect = bulk_baseline$effect,
    bootstrap_ci_lower = bulk_baseline_ci[1L],
    bootstrap_ci_upper = bulk_baseline_ci[2L],
    exact_p_two_sided = bulk_baseline$exact_p_two_sided,
    n_exact_assignments = bulk_baseline$n_assignments
  ),
  data.table(
    dataset = "GSE294385 paired macro-micro tumour",
    endpoint = "macro_minus_micro",
    n_units = nrow(spatial_rec_wide),
    n_rec_genes_defined = length(rec_genes),
    n_rec_genes_measured = spatial_rec$n_measured,
    rec_effect = spatial_baseline$effect,
    bootstrap_ci_lower = spatial_baseline_ci[1L],
    bootstrap_ci_upper = spatial_baseline_ci[2L],
    exact_p_two_sided = spatial_baseline$exact_p_two_sided,
    n_exact_assignments = spatial_baseline$n_assignments
  )
)

model_rows <- list()
score_rows <- list()
row_index <- 0L

for (rival_id in rival_registry$set_id) {
  rival_genes <- unique(gene_sets[[rival_id]])
  retained_rec_genes <- setdiff(rec_genes, rival_genes)

  bulk_rival <- score_program(bulk_gene_z, rival_genes)
  bulk_deoverlap <- score_program(bulk_gene_z, retained_rec_genes)
  spatial_rival <- score_program(spatial_gene_z, rival_genes)
  spatial_deoverlap <- score_program(spatial_gene_z, retained_rec_genes)
  if (any(vapply(
    list(bulk_rival, bulk_deoverlap, spatial_rival, spatial_deoverlap),
    function(x) is.null(x$score), logical(1)
  ))) {
    next
  }

  # Bulk endpoint.
  bulk_rival_score <- as.numeric(bulk_rival$score[bulk_metadata$sample_id])
  bulk_deoverlap_score <- as.numeric(bulk_deoverlap$score[bulk_metadata$sample_id])
  bulk_adjusted <- exact_adjusted_unpaired(
    bulk_rec_score, bulk_rival_score, bulk_group
  )
  bulk_deoverlap_unadjusted <- exact_two_group(bulk_deoverlap_score, bulk_group)
  bulk_deoverlap_adjusted <- exact_adjusted_unpaired(
    bulk_deoverlap_score, bulk_rival_score, bulk_group
  )
  bulk_rival_endpoint <- exact_two_group(bulk_rival_score, bulk_group)
  bulk_adjusted_seed <- stable_model_seed(
    "rHGP_minus_dHGP", rival_id, "adjusted_bootstrap"
  )
  set.seed(bulk_adjusted_seed)
  bulk_adjusted_ci <- bootstrap_unpaired_ci(
    bulk_rec_score, bulk_rival_score, bulk_group, bootstrap_iterations
  )
  bulk_deoverlap_seed <- stable_model_seed(
    "rHGP_minus_dHGP", rival_id, "deoverlap_bootstrap"
  )
  set.seed(bulk_deoverlap_seed)
  bulk_deoverlap_adjusted_ci <- bootstrap_unpaired_ci(
    bulk_deoverlap_score, bulk_rival_score, bulk_group, bootstrap_iterations
  )
  row_index <- row_index + 1L
  model_rows[[row_index]] <- data.table(
    dataset = "GSE151165 HGP-labelled bulk tumour",
    endpoint = "rHGP_minus_dHGP",
    rival_set_id = rival_id,
    n_units = nrow(bulk_metadata),
    n_rival_genes_defined = bulk_rival$n_declared,
    n_rival_genes_measured = bulk_rival$n_measured,
    n_rec_rival_overlap = length(intersect(rec_genes, rival_genes)),
    n_rec_genes_retained = length(retained_rec_genes),
    n_rec_genes_retained_measured = bulk_deoverlap$n_measured,
    rec_rival_score_correlation = cor(bulk_rec_score, bulk_rival_score),
    unadjusted_rec_effect = bulk_baseline$effect,
    adjusted_rec_effect = bulk_adjusted$effect,
    adjusted_bootstrap_ci_lower = bulk_adjusted_ci[1L],
    adjusted_bootstrap_ci_upper = bulk_adjusted_ci[2L],
    adjusted_bootstrap_seed = bulk_adjusted_seed,
    adjusted_exact_p_two_sided = bulk_adjusted$exact_p_two_sided,
    adjusted_effect_retention = bulk_adjusted$effect / bulk_baseline$effect,
    deoverlap_unadjusted_rec_effect = bulk_deoverlap_unadjusted$effect,
    deoverlap_unadjusted_effect_retention =
      bulk_deoverlap_unadjusted$effect / bulk_baseline$effect,
    deoverlap_adjusted_rec_effect = bulk_deoverlap_adjusted$effect,
    deoverlap_adjusted_bootstrap_ci_lower = bulk_deoverlap_adjusted_ci[1L],
    deoverlap_adjusted_bootstrap_ci_upper = bulk_deoverlap_adjusted_ci[2L],
    deoverlap_adjusted_bootstrap_seed = bulk_deoverlap_seed,
    deoverlap_adjusted_exact_p_two_sided =
      bulk_deoverlap_adjusted$exact_p_two_sided,
    deoverlap_adjusted_effect_retention =
      bulk_deoverlap_adjusted$effect / bulk_baseline$effect,
    rival_endpoint_effect = bulk_rival_endpoint$effect,
    rival_endpoint_exact_p_two_sided = bulk_rival_endpoint$exact_p_two_sided,
    n_exact_assignments = bulk_adjusted$n_assignments
  )

  score_rows[[paste0(rival_id, "__bulk")]] <- data.table(
    dataset = "GSE151165 HGP-labelled bulk tumour",
    unit = bulk_metadata$sample_id,
    group = bulk_metadata$hgp,
    rival_set_id = rival_id,
    rec_score = bulk_rec_score,
    rival_score = bulk_rival_score,
    deoverlap_rec_score = bulk_deoverlap_score
  )

  # Paired lesion-outgrowth endpoint.
  spatial_rival_dt <- copy(spatial_metadata)
  spatial_rival_dt[, `:=`(
    rival_score = as.numeric(spatial_rival$score[sample_id]),
    deoverlap_rec_score = as.numeric(spatial_deoverlap$score[sample_id])
  )]
  spatial_wide_scores <- dcast(
    spatial_rival_dt,
    patient ~ region,
    value.var = c("rival_score", "deoverlap_rec_score")
  )
  spatial_wide_scores <- spatial_rec_wide[
    spatial_wide_scores,
    on = "patient"
  ]
  spatial_wide_scores[, `:=`(
    rival_difference = rival_score_macro_tumour - rival_score_micro_tumour,
    deoverlap_rec_difference =
      deoverlap_rec_score_macro_tumour - deoverlap_rec_score_micro_tumour
  )]
  spatial_adjusted <- exact_adjusted_paired(
    spatial_wide_scores$rec_difference,
    spatial_wide_scores$rival_difference
  )
  spatial_deoverlap_unadjusted <- exact_paired(
    spatial_wide_scores$deoverlap_rec_difference
  )
  spatial_deoverlap_adjusted <- exact_adjusted_paired(
    spatial_wide_scores$deoverlap_rec_difference,
    spatial_wide_scores$rival_difference
  )
  spatial_rival_endpoint <- exact_paired(spatial_wide_scores$rival_difference)
  spatial_adjusted_seed <- stable_model_seed(
    "macro_minus_micro", rival_id, "adjusted_bootstrap"
  )
  set.seed(spatial_adjusted_seed)
  spatial_adjusted_ci <- bootstrap_paired_ci(
    spatial_wide_scores$rec_difference,
    spatial_wide_scores$rival_difference,
    bootstrap_iterations
  )
  spatial_deoverlap_seed <- stable_model_seed(
    "macro_minus_micro", rival_id, "deoverlap_bootstrap"
  )
  set.seed(spatial_deoverlap_seed)
  spatial_deoverlap_adjusted_ci <- bootstrap_paired_ci(
    spatial_wide_scores$deoverlap_rec_difference,
    spatial_wide_scores$rival_difference,
    bootstrap_iterations
  )
  row_index <- row_index + 1L
  model_rows[[row_index]] <- data.table(
    dataset = "GSE294385 paired macro-micro tumour",
    endpoint = "macro_minus_micro",
    rival_set_id = rival_id,
    n_units = nrow(spatial_wide_scores),
    n_rival_genes_defined = spatial_rival$n_declared,
    n_rival_genes_measured = spatial_rival$n_measured,
    n_rec_rival_overlap = length(intersect(rec_genes, rival_genes)),
    n_rec_genes_retained = length(retained_rec_genes),
    n_rec_genes_retained_measured = spatial_deoverlap$n_measured,
    rec_rival_score_correlation = cor(
      spatial_wide_scores$rec_difference,
      spatial_wide_scores$rival_difference
    ),
    unadjusted_rec_effect = spatial_baseline$effect,
    adjusted_rec_effect = spatial_adjusted$effect,
    adjusted_bootstrap_ci_lower = spatial_adjusted_ci[1L],
    adjusted_bootstrap_ci_upper = spatial_adjusted_ci[2L],
    adjusted_bootstrap_seed = spatial_adjusted_seed,
    adjusted_exact_p_two_sided = spatial_adjusted$exact_p_two_sided,
    adjusted_effect_retention = spatial_adjusted$effect / spatial_baseline$effect,
    deoverlap_unadjusted_rec_effect = spatial_deoverlap_unadjusted$effect,
    deoverlap_unadjusted_effect_retention =
      spatial_deoverlap_unadjusted$effect / spatial_baseline$effect,
    deoverlap_adjusted_rec_effect = spatial_deoverlap_adjusted$effect,
    deoverlap_adjusted_bootstrap_ci_lower = spatial_deoverlap_adjusted_ci[1L],
    deoverlap_adjusted_bootstrap_ci_upper = spatial_deoverlap_adjusted_ci[2L],
    deoverlap_adjusted_bootstrap_seed = spatial_deoverlap_seed,
    deoverlap_adjusted_exact_p_two_sided =
      spatial_deoverlap_adjusted$exact_p_two_sided,
    deoverlap_adjusted_effect_retention =
      spatial_deoverlap_adjusted$effect / spatial_baseline$effect,
    rival_endpoint_effect = spatial_rival_endpoint$effect,
    rival_endpoint_exact_p_two_sided =
      spatial_rival_endpoint$exact_p_two_sided,
    n_exact_assignments = spatial_adjusted$n_assignments
  )

  score_rows[[paste0(rival_id, "__spatial")]] <- data.table(
    dataset = "GSE294385 paired macro-micro tumour",
    unit = spatial_wide_scores$patient,
    group = "macro_minus_micro",
    rival_set_id = rival_id,
    rec_score = spatial_wide_scores$rec_difference,
    rival_score = spatial_wide_scores$rival_difference,
    deoverlap_rec_score = spatial_wide_scores$deoverlap_rec_difference
  )
}

models <- rbindlist(model_rows, use.names = TRUE, fill = TRUE)
models <- rival_registry[
  models,
  on = c(set_id = "rival_set_id")
]
setnames(models, "set_id", "rival_set_id")
models[, `:=`(
  adjusted_fdr_bh = p.adjust(adjusted_exact_p_two_sided, method = "BH"),
  deoverlap_adjusted_fdr_bh = p.adjust(
    deoverlap_adjusted_exact_p_two_sided,
    method = "BH"
  ),
  rival_endpoint_fdr_bh = p.adjust(
    rival_endpoint_exact_p_two_sided,
    method = "BH"
  )
), by = dataset]
setorder(models, dataset, role, rival_set_id)
scores <- rbindlist(score_rows, use.names = TRUE, fill = TRUE)

fwrite(
  baseline,
  file.path(out_dir, "rec_external_rival_baseline.tsv"),
  sep = "\t"
)
fwrite(
  models,
  file.path(out_dir, "rec_external_rival_models.tsv"),
  sep = "\t"
)
fwrite(
  scores,
  file.path(out_dir, "rec_external_rival_scores.tsv.gz"),
  sep = "\t"
)
capture.output(
  sessionInfo(),
  file = file.path(out_dir, "rec_external_rival_sessionInfo.txt")
)

print(baseline)
focus_ids <- c(
  "OGDEN_HYPOXIA_MP6", "HALLMARK_HYPOXIA",
  "HALLMARK_UNFOLDED_PROTEIN_RESPONSE", "HALLMARK_E2F_TARGETS",
  "HALLMARK_G2M_CHECKPOINT", "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_FATTY_ACID_METABOLISM", "HALLMARK_GLYCOLYSIS",
  "PENGWINKLER_PALMITATE_UP_TOP100", "PENGWINKLER_CAUSAL_AXIS_5",
  "REACTOME_COLLAGEN_BIOSYNTHESIS_AND_MODIFYING_ENZYMES", "PAN_EPITHELIAL_7",
  "HEPATOCYTE_CONTEXT_6"
)
print(models[
  rival_set_id %in% focus_ids,
  .(dataset, rival_set_id, role, n_rec_rival_overlap,
    rec_rival_score_correlation, adjusted_rec_effect,
    adjusted_bootstrap_ci_lower, adjusted_bootstrap_ci_upper,
    adjusted_exact_p_two_sided, adjusted_fdr_bh,
    deoverlap_adjusted_rec_effect, deoverlap_adjusted_effect_retention,
    rival_endpoint_effect, rival_endpoint_exact_p_two_sided)
])
