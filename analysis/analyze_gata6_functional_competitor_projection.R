#!/usr/bin/env Rscript

# Functional competitor projection: GATA6 loss versus frozen REC-state programmes
# Date: 2026-09-01
# Random seed: 42
# Specification: metadata/gata6_functional_competitor_projection_spec_2026-09-01.md

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
})

options(stringsAsFactors = FALSE)
random_seed <- 42L
bootstrap_iterations <- 10000L
set.seed(random_seed)

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
root <- if (length(file_arg)) {
  normalizePath(file.path(dirname(sub("^--file=", "", file_arg[[1L]])), ".."), mustWork = TRUE)
} else {
  normalizePath(".", mustWork = TRUE)
}

phase2_dir <- file.path(
  root, "analysis_results", "deep_biology_upgrade_2026-08-31",
  "phase2_mechanistic_specificity"
)
phase3_dir <- file.path(
  root, "analysis_results", "deep_biology_upgrade_2026-08-31",
  "phase3_external_regulatory_projection"
)
out_dir <- file.path(
  root, "analysis_results", "deep_biology_upgrade_2026-08-31",
  "phase5_gata6_functional_competitor_projection"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

counts_path <- file.path(
  root, "data_sources", "Goto_2026_GATA6", "GSE290753",
  "GSE290753_ngoto_RNAseq_GATA6KO_Ct.txt.gz"
)
series_path <- file.path(
  root, "data_sources", "Goto_2026_GATA6", "GSE290753",
  "GSE290753_series_matrix.txt.gz"
)
gtf_path <- file.path(
  root, "data_sources", "Goto_2026_GATA6", "reference",
  "Mus_musculus.GRCm38.101.gtf.gz"
)
phase2_sets_path <- file.path(phase2_dir, "frozen_gene_sets_long.tsv")
phase3_sets_path <- file.path(
  phase3_dir, "frozen_external_regulatory_gene_sets_long.tsv"
)
coexistence_path <- file.path(
  phase2_dir, "ogden_rec_coexistence_locked_genes.tsv"
)
spec_path <- file.path(
  root, "metadata", "gata6_functional_competitor_projection_spec_2026-09-01.md"
)
required <- c(
  counts_path, series_path, gtf_path, phase2_sets_path, phase3_sets_path,
  coexistence_path, spec_path
)
if (any(!file.exists(required))) stop("One or more GATA6-projection inputs are missing")

tol <- sqrt(.Machine$double.eps)

write_tsv <- function(x, filename) {
  fwrite(as.data.table(x), file.path(out_dir, filename), sep = "\t", na = "NA")
}

parse_gtf_gene_map <- function(path) {
  command <- sprintf(
    "gzip -dc %s | awk -F '\\t' '$3 == \"gene\" {print $9}'",
    shQuote(path)
  )
  attributes <- fread(cmd = command, header = FALSE, sep = "\n", quote = "")$V1
  gene_id <- sub('.*gene_id "([^"]+)".*', "\\1", attributes)
  gene_name <- sub('.*gene_name "([^"]+)".*', "\\1", attributes)
  result <- data.table(
    ensembl_gene_id = sub("\\.[0-9]+$", "", gene_id),
    gene = toupper(gene_name)
  )
  result <- result[
    !is.na(ensembl_gene_id) & ensembl_gene_id != "" &
      !is.na(gene) & gene != ""
  ]
  unique(result, by = "ensembl_gene_id")
}

safe_gene_z <- function(x) {
  means <- rowMeans(x, na.rm = TRUE)
  sds <- apply(x, 1L, sd, na.rm = TRUE)
  keep <- is.finite(sds) & sds > 0
  z <- sweep(x[keep, , drop = FALSE], 1L, means[keep], "-")
  sweep(z, 1L, sds[keep], "/")
}

hedges_g <- function(values, group) {
  knockout <- values[group == 1L]
  control <- values[group == 0L]
  n_knockout <- length(knockout)
  n_control <- length(control)
  pooled_variance <- (
    (n_knockout - 1) * var(knockout) +
      (n_control - 1) * var(control)
  ) / (n_knockout + n_control - 2)
  if (!is.finite(pooled_variance) || pooled_variance <= 0) return(NA_real_)
  correction <- 1 - 3 / (4 * (n_knockout + n_control) - 9)
  correction * (mean(knockout) - mean(control)) / sqrt(pooled_variance)
}

exact_two_group <- function(values, group) {
  n_knockout <- sum(group == 1L)
  observed <- mean(values[group == 1L]) - mean(values[group == 0L])
  allocations <- combn(seq_along(values), n_knockout)
  knockout_sums <- colSums(matrix(values[allocations], nrow = n_knockout))
  null <- knockout_sums / n_knockout -
    (sum(values) - knockout_sums) / (length(values) - n_knockout)
  list(
    effect = observed,
    p_two_sided = mean(abs(null) >= abs(observed) - tol),
    n_allocations = length(null)
  )
}

bootstrap_group_difference <- function(values, group, iterations) {
  knockout_index <- which(group == 1L)
  control_index <- which(group == 0L)
  estimates <- replicate(iterations, {
    mean(values[sample(knockout_index, length(knockout_index), replace = TRUE)]) -
      mean(values[sample(control_index, length(control_index), replace = TRUE)])
  })
  as.numeric(quantile(estimates, c(0.025, 0.975), names = FALSE, type = 7))
}

phase2_sets <- fread(phase2_sets_path)
phase3_sets <- fread(phase3_sets_path)
coexistence <- fread(coexistence_path)
phase2_sets[, gene := toupper(gene)]
phase3_sets[, gene := toupper(gene)]
coexistence[, gene := toupper(gene)]

# data.table's scoping makes a short explicit helper clearer than non-standard evaluation.
get_set <- function(table, identifier) {
  unique(as.character(table[table$set_id == identifier, "gene"][[1L]]))
}

gene_sets <- list(
  REC_PROGRAM_TOP50 = get_set(phase2_sets, "FROZEN_REC_TOP50"),
  CORE_HRC = get_set(phase2_sets, "CANELLAS_CORE_HRC"),
  PUBLISHED_FETAL = get_set(phase2_sets, "PUBLISHED_FETAL"),
  PARTIAL_EMT = get_set(phase2_sets, "PUBLISHED_PEMT"),
  HALLMARK_EMT = get_set(
    phase2_sets, "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"
  ),
  REACTOME_TIGHT_JUNCTION_INTERACTIONS = get_set(
    phase3_sets, "REACTOME_TIGHT_JUNCTION_INTERACTIONS"
  ),
  JUNB_TIGHT_JUNCTION_TARGETS_5 = get_set(
    phase3_sets, "JUNB_TIGHT_JUNCTION_TARGETS_5"
  ),
  ACTIN_TURNOVER = unique(as.character(coexistence[
    program_id == "ACTIN_TURNOVER_UNION" & gene_role == "target", gene
  ])),
  HALLMARK_E2F_TARGETS = get_set(phase2_sets, "HALLMARK_E2F_TARGETS"),
  HALLMARK_G2M_CHECKPOINT = get_set(phase2_sets, "HALLMARK_G2M_CHECKPOINT")
)
if (any(lengths(gene_sets) == 0L)) stop("At least one frozen programme is empty")

named_genes <- c(
  "GATA6", "HNF4A", "LGR5", "MKI67", "EMP1", "KLF4", "EPCAM", "ZEB1",
  "CDH17", "TJP1", "CLDN2", "CLDN3", "CLDN4", "CLDN7", "CRB3", "F11R"
)

counts_command <- sprintf("gzip -dc %s", shQuote(counts_path))
counts_raw <- fread(cmd = counts_command)
setnames(counts_raw, 1L, "ensembl_gene_id")
counts_raw[, ensembl_gene_id := sub("\\.[0-9]+$", "", ensembl_gene_id)]
sample_ids <- setdiff(names(counts_raw), "ensembl_gene_id")
expected_samples <- c(
  "Control_replicate1", "Control_replicate2", "Control_replicate3",
  "GATA6KO_replicate1", "GATA6KO_replicate2", "GATA6KO_replicate3"
)
if (!identical(sample_ids, expected_samples)) stop("Unexpected GSE290753 sample order")
if (any(vapply(counts_raw[, ..sample_ids], function(x) any(x < 0), logical(1L)))) {
  stop("Negative values found in the declared count matrix")
}
if (any(vapply(counts_raw[, ..sample_ids], function(x) any(x != floor(x)), logical(1L)))) {
  stop("GSE290753 processed matrix is not integer count data")
}

gene_map <- parse_gtf_gene_map(gtf_path)
counts_annotated <- merge(counts_raw, gene_map, by = "ensembl_gene_id", all.x = TRUE)
mapping_qc <- data.table(
  n_input_ensembl_ids = nrow(counts_raw),
  n_mapped_ensembl_ids = sum(!is.na(counts_annotated$gene)),
  n_unmapped_ensembl_ids = sum(is.na(counts_annotated$gene)),
  n_unique_mapped_symbols = uniqueN(counts_annotated$gene, na.rm = TRUE)
)
counts_annotated <- counts_annotated[!is.na(gene) & gene != ""]
counts_by_gene <- counts_annotated[
  , lapply(.SD, sum), by = gene, .SDcols = sample_ids
]
count_matrix <- as.matrix(counts_by_gene[, ..sample_ids])
storage.mode(count_matrix) <- "integer"
rownames(count_matrix) <- counts_by_gene$gene

condition <- factor(
  c(rep("control", 3L), rep("GATA6_KO", 3L)),
  levels = c("control", "GATA6_KO")
)
group <- as.integer(condition == "GATA6_KO")
sample_metadata <- data.table(
  sample_id = sample_ids,
  condition = as.character(condition),
  raw_library_size = colSums(count_matrix),
  detected_genes = colSums(count_matrix > 0)
)

dge <- DGEList(count_matrix, group = condition)
keep <- filterByExpr(dge, group = condition)
dge <- calcNormFactors(dge[keep, , keep.lib.sizes = FALSE])
design <- model.matrix(~condition)
voom_fit <- voom(dge, design, plot = FALSE)
fit <- eBayes(lmFit(voom_fit, design))
log_cpm <- voom_fit$E
gene_z <- safe_gene_z(log_cpm)

module_score_rows <- list()
module_coverage_rows <- list()
module_contrast_rows <- list()
loo_rows <- list()
module_scores <- list()
for (set_id in names(gene_sets)) {
  declared <- unique(toupper(gene_sets[[set_id]]))
  measured <- intersect(declared, rownames(gene_z))
  if (!length(measured)) stop(paste("No measurable genes for", set_id))
  values <- colMeans(gene_z[measured, , drop = FALSE])
  module_scores[[set_id]] <- values
  module_score_rows[[set_id]] <- data.table(
    sample_id = sample_ids,
    condition = as.character(condition),
    set_id = set_id,
    score = as.numeric(values[sample_ids])
  )
  module_coverage_rows[[set_id]] <- data.table(
    set_id = set_id,
    n_declared = length(declared),
    n_measured_nonconstant = length(measured),
    coverage_fraction = length(measured) / length(declared),
    measured_genes = paste(measured, collapse = ";"),
    missing_or_filtered_genes = paste(setdiff(declared, measured), collapse = ";")
  )
  values <- as.numeric(values[sample_ids])
  exact <- exact_two_group(values, group)
  interval <- bootstrap_group_difference(values, group, bootstrap_iterations)
  module_contrast_rows[[set_id]] <- data.table(
    set_id = set_id,
    n_libraries = length(values),
    n_gata6_ko = sum(group == 1L),
    n_control = sum(group == 0L),
    gata6_ko_mean = mean(values[group == 1L]),
    control_mean = mean(values[group == 0L]),
    gata6_ko_minus_control = exact$effect,
    bootstrap_ci_lower = interval[[1L]],
    bootstrap_ci_upper = interval[[2L]],
    hedges_g = hedges_g(values, group),
    ko_above_control_median_n = sum(values[group == 1L] > median(values[group == 0L])),
    exact_p_two_sided = exact$p_two_sided,
    n_exact_allocations = exact$n_allocations
  )
  if (set_id %in% c(
    "REC_PROGRAM_TOP50", "REACTOME_TIGHT_JUNCTION_INTERACTIONS",
    "JUNB_TIGHT_JUNCTION_TARGETS_5"
  )) {
    loo_rows[[set_id]] <- rbindlist(lapply(measured, function(omitted_gene) {
      retained <- setdiff(measured, omitted_gene)
      retained_values <- colMeans(gene_z[retained, , drop = FALSE])
      data.table(
        set_id = set_id,
        omitted_gene = omitted_gene,
        retained_gene_count = length(retained),
        gata6_ko_minus_control =
          mean(retained_values[group == 1L]) - mean(retained_values[group == 0L])
      )
    }))
  }
}

module_scores_table <- rbindlist(module_score_rows, use.names = TRUE)
module_coverage <- rbindlist(module_coverage_rows, use.names = TRUE)
module_contrasts <- rbindlist(module_contrast_rows, use.names = TRUE)
module_contrasts[, fdr_bh_module_family := p.adjust(exact_p_two_sided, "BH")]
module_loo <- rbindlist(loo_rows, use.names = TRUE)

camera_index <- lapply(gene_sets, function(genes) {
  which(rownames(voom_fit$E) %chin% unique(toupper(genes)))
})
camera_result <- as.data.table(camera(
  voom_fit,
  index = camera_index,
  design = design,
  contrast = ncol(design),
  sort = FALSE
), keep.rownames = "set_id")
setnames(
  camera_result,
  c("NGenes", "Direction", "PValue", "FDR"),
  c("n_genes", "direction", "p_value", "fdr_bh_camera_family")
)

gene_table <- as.data.table(topTable(
  fit, coef = ncol(design), number = Inf, sort.by = "none"
), keep.rownames = "gene")
gene_table[, gene := toupper(gene)]
moderated_se <- fit$stdev.unscaled[, ncol(design)] * sqrt(fit$s2.post)
critical_t <- qt(0.975, df = fit$df.total)
ci_table <- data.table(
  gene = toupper(rownames(fit$coefficients)),
  ci_lower = fit$coefficients[, ncol(design)] - critical_t * moderated_se,
  ci_upper = fit$coefficients[, ncol(design)] + critical_t * moderated_se
)
named_effects <- merge(
  data.table(gene_order = seq_along(named_genes), gene = named_genes),
  gene_table[, .(
    gene,
    log2fc_gata6_ko_minus_control = logFC,
    average_log2_expression = AveExpr,
    moderated_t = t,
    p_value = P.Value,
    fdr_bh_all_tested_genes = adj.P.Val
  )],
  by = "gene", all.x = TRUE, sort = FALSE
)
named_effects <- merge(named_effects, ci_table, by = "gene", all.x = TRUE, sort = FALSE)
setorder(named_effects, gene_order)
named_effects[, fdr_bh_named_family := p.adjust(p_value, "BH")]

sample_metadata[, `:=`(
  tmm_norm_factor = dge$samples$norm.factors,
  effective_library_size = dge$samples$lib.size * dge$samples$norm.factors
)]
sample_correlation <- as.data.table(cor(log_cpm, method = "pearson"), keep.rownames = "sample_id")

write_tsv(mapping_qc, "gse290753_gene_mapping_qc.tsv")
write_tsv(sample_metadata, "gse290753_sample_qc.tsv")
write_tsv(sample_correlation, "gse290753_sample_correlation.tsv")
write_tsv(module_coverage, "gse290753_programme_coverage.tsv")
write_tsv(module_scores_table, "gse290753_programme_scores.tsv")
write_tsv(module_contrasts, "gse290753_programme_contrasts.tsv")
write_tsv(module_loo, "gse290753_programme_leave_one_gene_out.tsv")
write_tsv(camera_result, "gse290753_camera_results.tsv")
write_tsv(named_effects, "gse290753_named_gene_effects.tsv")
fwrite(
  gene_table,
  file.path(out_dir, "gse290753_limma_all_genes.tsv.gz"),
  sep = "\t", na = "NA", compress = "gzip"
)

summary_lines <- c(
  "# GSE290753 GATA6-loss functional competitor projection",
  "",
  "Three GATA6-knockout and three control organoid RNA-seq libraries were analysed as independent libraries. Exact label-allocation P values have only 20 possible allocations and are therefore coarse.",
  "",
  "## Frozen programme contrasts",
  ""
)
for (current_set in names(gene_sets)) {
  row <- module_contrasts[module_contrasts$set_id == current_set]
  camera_row <- camera_result[camera_result$set_id == current_set]
  summary_lines <- c(
    summary_lines,
    sprintf(
      paste0(
        "- `%s`: KO-control %+.3f (95%% library-bootstrap interval %+.3f to %+.3f; ",
        "Hedges' g %+.3f; exact P=%.3g); camera %s, P=%.3g, FDR=%.3g."
      ),
      current_set, row$gata6_ko_minus_control, row$bootstrap_ci_lower,
      row$bootstrap_ci_upper, row$hedges_g, row$exact_p_two_sided,
      camera_row$direction, camera_row$p_value,
      camera_row$fdr_bh_camera_family
    )
  )
}
summary_lines <- c(
  summary_lines,
  "",
  "## Result interpretation",
  "",
  "The perturbation behaves as expected at its defining genes: GATA6, HNF4A and LGR5 all decrease. Despite that LGR5-low direction, the frozen REC and core-HRC programmes, broad tight-junction programme and focused five-gene junction subset all decrease. Every leave-one-gene-out REC and junction estimate remains negative, and E2F/G2M are not reduced. GATA6 loss is therefore a strong functional counterexample: source-described metastatic lineage plasticity and LGR5 loss are not sufficient to produce the cohesive low-cycle REC phenotype.",
  "",
  "The cross-study fetal and partial-EMT sets used here are not the source paper's own GATA6-response signatures. Their directions should not be used to dispute the source study's fetal/basal conclusions; their role is only to test transfer of the already frozen manuscript programmes.",
  "",
  "## Interpretation boundary",
  "",
  "This is a functional competitor projection without HGP labels. Concordance can show overlap with a GATA6-loss metastatic-plasticity route; discordance can distinguish the replacement-associated state from generic LGR5-low plasticity. Neither result establishes GATA6 as an rHGP driver.",
  ""
)
writeLines(summary_lines, file.path(out_dir, "analysis_summary.md"), useBytes = TRUE)

run_info <- data.table(
  analysis_date = "2026-09-01",
  random_seed = random_seed,
  bootstrap_iterations = bootstrap_iterations,
  r_version = R.version.string,
  edgeR_version = as.character(packageVersion("edgeR")),
  limma_version = as.character(packageVersion("limma")),
  inferential_unit = "RNA-seq library",
  n_control = sum(group == 0L),
  n_gata6_ko = sum(group == 1L),
  exact_allocations = choose(length(group), sum(group == 1L)),
  source_annotation = "Ensembl release 101, mm10/GRCm38",
  counts_md5 = unname(tools::md5sum(counts_path)),
  series_md5 = unname(tools::md5sum(series_path)),
  gtf_md5 = unname(tools::md5sum(gtf_path)),
  phase2_gene_sets_md5 = unname(tools::md5sum(phase2_sets_path)),
  phase3_gene_sets_md5 = unname(tools::md5sum(phase3_sets_path)),
  specification_md5 = unname(tools::md5sum(spec_path))
)
write_tsv(run_info, "gse290753_run_info.tsv")
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))

manifest <- data.table(
  file = sort(list.files(out_dir, full.names = FALSE)),
  role = "analysis_output"
)
write_tsv(manifest, "output_manifest.tsv")

cat("Wrote GATA6 functional competitor projection to", out_dir, "\n")
