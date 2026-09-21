#!/usr/bin/env Rscript

# Analysis: Paired Ogden REC-versus-epithelial-state pseudobulk contrasts
# Date: 2026-09-01
# Random seed: 42
# R and package versions are written to sessionInfo output.

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
})

set.seed(42)

root <- normalizePath(getwd(), mustWork = TRUE)
out_dir <- file.path(
  root,
  "analysis_results",
  "deep_biology_upgrade_2026-08-31",
  "phase2_mechanistic_specificity"
)
counts_path <- file.path(out_dir, "ogden_direct_state_pseudobulk_counts.tsv.gz")
metadata_path <- file.path(out_dir, "ogden_direct_state_pseudobulk_metadata.tsv")
eligibility_path <- file.path(out_dir, "ogden_direct_contrast_eligibility.tsv")
gene_set_path <- file.path(out_dir, "frozen_gene_sets_long.tsv")
registry_path <- file.path(out_dir, "frozen_gene_set_registry.tsv")

required_paths <- c(
  counts_path,
  metadata_path,
  eligibility_path,
  gene_set_path,
  registry_path
)
if (any(!file.exists(required_paths))) {
  stop("One or more required Phase 2 inputs are missing")
}

exact_signflip_p <- function(differences) {
  differences <- as.numeric(differences)
  n <- length(differences)
  if (n < 1L || n > 20L || any(!is.finite(differences))) {
    return(NA_real_)
  }
  observed <- abs(mean(differences))
  patterns <- 0:(2^n - 1L)
  null_values <- vapply(
    patterns,
    function(pattern) {
      signs <- ifelse(bitwAnd(pattern, bitwShiftL(1L, 0:(n - 1L))) > 0L, 1, -1)
      abs(mean(differences * signs))
    },
    numeric(1)
  )
  mean(null_values >= observed - sqrt(.Machine$double.eps))
}

bootstrap_mean_ci <- function(values, iterations = 10000L) {
  values <- as.numeric(values)
  n <- length(values)
  if (n < 2L || any(!is.finite(values))) {
    return(c(lower = NA_real_, upper = NA_real_))
  }
  estimates <- replicate(
    iterations,
    mean(sample(values, size = n, replace = TRUE))
  )
  quantile(estimates, probs = c(0.025, 0.975), names = FALSE, type = 6)
}

counts_dt <- fread(cmd = paste("gzip -dc", shQuote(counts_path)))
if (!"gene" %in% names(counts_dt) || anyDuplicated(counts_dt$gene)) {
  stop("Pseudobulk count table requires unique gene symbols")
}
genes <- toupper(as.character(counts_dt$gene))
count_matrix <- as.matrix(counts_dt[, -"gene"])
storage.mode(count_matrix) <- "double"
rownames(count_matrix) <- genes
if (any(!is.finite(count_matrix)) || any(count_matrix < 0)) {
  stop("Counts must be finite and non-negative")
}

metadata <- fread(metadata_path)
eligibility <- fread(eligibility_path)
gene_sets_long <- fread(gene_set_path)
registry <- fread(registry_path)
if (!all(metadata$sample_id %in% colnames(count_matrix))) {
  stop("Metadata sample IDs are not all present in the count matrix")
}

eligible_comparators <- eligibility[
  eligible_minimum_five_patients == TRUE,
  comparator_state
]
if (length(eligible_comparators) == 0L) {
  stop("No predeclared direct contrast is eligible")
}

all_gene_results <- list()
contrast_summaries <- list()
program_values <- list()
program_effects <- list()

for (comparator in eligible_comparators) {
  paired_patients <- intersect(
    metadata[state == "REC" & eligible_min_cells == TRUE, patient],
    metadata[state == comparator & eligible_min_cells == TRUE, patient]
  )
  paired_patients <- sort(unique(paired_patients))
  selected <- metadata[
    patient %in% paired_patients & state %in% c(comparator, "REC") &
      eligible_min_cells == TRUE
  ]
  setorder(selected, patient, state)
  if (nrow(selected) != 2L * length(paired_patients) ||
      any(selected[, .N, by = patient]$N != 2L)) {
    stop(sprintf("Incomplete pairing for REC versus %s", comparator))
  }

  selected_counts <- count_matrix[, selected$sample_id, drop = FALSE]
  selected$patient <- factor(selected$patient, levels = paired_patients)
  selected$state <- factor(selected$state, levels = c(comparator, "REC"))
  design <- model.matrix(~ patient + state, data = selected)
  if (qr(design)$rank != ncol(design) || !"stateREC" %in% colnames(design)) {
    stop(sprintf("Invalid paired design for REC versus %s", comparator))
  }

  dge <- DGEList(counts = selected_counts)
  keep <- filterByExpr(dge, design = design)
  if (sum(keep) < 1000L) {
    stop(sprintf("Too few expressed genes retained for REC versus %s", comparator))
  }
  dge <- calcNormFactors(dge[keep, , keep.lib.sizes = FALSE], method = "TMM")
  voom_fit <- voom(dge, design = design, plot = FALSE)
  fit <- eBayes(lmFit(voom_fit, design), robust = TRUE)
  table <- as.data.table(topTable(
    fit,
    coef = "stateREC",
    number = Inf,
    sort.by = "none",
    adjust.method = "BH"
  ), keep.rownames = "gene")
  setnames(
    table,
    c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B"),
    c("log2fc_rec_minus_comparator", "average_log2_expression", "moderated_t",
      "p_value", "fdr_bh", "log_odds_de")
  )
  table[, `:=`(
    comparator_state = comparator,
    n_paired_patients = length(paired_patients),
    contrast = paste0("REC_vs_", gsub(" ", "_", comparator))
  )]
  setcolorder(
    table,
    c("contrast", "comparator_state", "n_paired_patients", "gene",
      "log2fc_rec_minus_comparator", "average_log2_expression", "moderated_t",
      "p_value", "fdr_bh", "log_odds_de")
  )
  all_gene_results[[comparator]] <- table

  contrast_summaries[[comparator]] <- data.table(
    contrast = paste0("REC_vs_", gsub(" ", "_", comparator)),
    comparator_state = comparator,
    n_paired_patients = length(paired_patients),
    paired_patient_ids = paste(paired_patients, collapse = ";"),
    n_rec_cells = selected[state == "REC", sum(n_cells)],
    n_comparator_cells = selected[state == comparator, sum(n_cells)],
    n_tested_genes = nrow(table),
    n_fdr05_rec_up = table[fdr_bh < 0.05 & log2fc_rec_minus_comparator > 0, .N],
    n_fdr05_rec_down = table[fdr_bh < 0.05 & log2fc_rec_minus_comparator < 0, .N]
  )

  log_cpm <- cpm(dge, log = TRUE, prior.count = 2)
  gene_sd <- apply(log_cpm, 1L, sd)
  valid_z <- is.finite(gene_sd) & gene_sd > 0
  z_matrix <- matrix(NA_real_, nrow = nrow(log_cpm), ncol = ncol(log_cpm),
                     dimnames = dimnames(log_cpm))
  z_matrix[valid_z, ] <- t(scale(t(log_cpm[valid_z, , drop = FALSE])))

  for (set_name in registry$set_id) {
    declared_genes <- gene_sets_long[set_id == set_name, unique(toupper(gene))]
    present_genes <- intersect(declared_genes, rownames(z_matrix)[valid_z])
    if (length(present_genes) < 3L) {
      next
    }
    scores <- colMeans(z_matrix[present_genes, , drop = FALSE], na.rm = TRUE)
    values <- data.table(
      sample_id = names(scores),
      score = as.numeric(scores)
    )[selected, on = "sample_id"]
    values[, `:=`(
      set_id = set_name,
      comparator_state = comparator,
      n_declared_genes = length(declared_genes),
      n_scored_genes = length(present_genes)
    )]
    program_values[[paste(comparator, set_name, sep = "__")]] <- values[, .(
      comparator_state, set_id, patient = as.character(patient), state = as.character(state),
      sample_id, n_cells, n_declared_genes, n_scored_genes, score
    )]

    paired_scores <- dcast(
      values,
      patient ~ state,
      value.var = "score"
    )
    paired_scores[, difference_rec_minus_comparator := REC - get(comparator)]
    differences <- paired_scores$difference_rec_minus_comparator
    ci <- bootstrap_mean_ci(differences)
    difference_sd <- sd(differences)
    program_effects[[paste(comparator, set_name, sep = "__")]] <- data.table(
      comparator_state = comparator,
      set_id = set_name,
      n_paired_patients = length(differences),
      n_declared_genes = length(declared_genes),
      n_scored_genes = length(present_genes),
      mean_difference = mean(differences),
      median_difference = median(differences),
      bootstrap_ci_lower = ci[1L],
      bootstrap_ci_upper = ci[2L],
      paired_standardized_effect = ifelse(
        is.finite(difference_sd) && difference_sd > 0,
        mean(differences) / difference_sd,
        NA_real_
      ),
      patient_fraction_positive = mean(differences > 0),
      exact_signflip_p = exact_signflip_p(differences)
    )
  }
}

all_gene_dt <- rbindlist(all_gene_results, use.names = TRUE)
contrast_dt <- rbindlist(contrast_summaries, use.names = TRUE)
program_value_dt <- rbindlist(program_values, use.names = TRUE)
program_effect_dt <- rbindlist(program_effects, use.names = TRUE)
program_effect_dt <- registry[
  program_effect_dt,
  on = "set_id"
]
program_effect_dt[, fdr_bh_within_contrast := p.adjust(exact_signflip_p, method = "BH"),
                  by = comparator_state]
setorder(program_effect_dt, comparator_state, exact_signflip_p, set_id)

fwrite(
  all_gene_dt,
  file.path(out_dir, "ogden_direct_state_limma_all.tsv.gz"),
  sep = "\t"
)
fwrite(
  contrast_dt,
  file.path(out_dir, "ogden_direct_state_contrast_summary.tsv"),
  sep = "\t"
)
fwrite(
  program_value_dt,
  file.path(out_dir, "ogden_direct_state_program_scores.tsv.gz"),
  sep = "\t"
)
fwrite(
  program_effect_dt,
  file.path(out_dir, "ogden_direct_state_program_effects.tsv"),
  sep = "\t"
)

capture.output(
  sessionInfo(),
  file = file.path(out_dir, "ogden_direct_state_sessionInfo.txt")
)

print(contrast_dt)
print(program_effect_dt[
  set_id %in% c(
    "FROZEN_REC_TOP50", "HALLMARK_HYPOXIA",
    "HALLMARK_UNFOLDED_PROTEIN_RESPONSE", "HALLMARK_E2F_TARGETS",
    "HALLMARK_MYC_TARGETS_V1", "HALLMARK_FATTY_ACID_METABOLISM",
    "PENGWINKLER_PALMITATE_UP_TOP100", "PENGWINKLER_CAUSAL_AXIS_5"
  ),
  .(comparator_state, set_id, n_paired_patients, n_scored_genes,
    mean_difference, bootstrap_ci_lower, bootstrap_ci_upper,
    patient_fraction_positive, exact_signflip_p, fdr_bh_within_contrast)
])
