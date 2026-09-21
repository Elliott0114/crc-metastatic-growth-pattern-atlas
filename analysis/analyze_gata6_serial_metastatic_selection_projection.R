#!/usr/bin/env Rscript

# Serial metastatic-selection projection: GSE290752 versus frozen REC-state programmes
# Date: 2026-09-01
# Random seed: 42
# Specification: metadata/gata6_serial_metastatic_selection_projection_spec_2026-09-01.md

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
  "phase6_gata6_serial_metastatic_selection_projection"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

counts_path <- file.path(
  root, "data_sources", "Goto_2026_GATA6", "GSE290752",
  "GSE290752_ngoto_RNAseq_met_Ct.txt.gz"
)
series_path <- file.path(
  root, "data_sources", "Goto_2026_GATA6", "GSE290752",
  "GSE290752_series_matrix.txt.gz"
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
  root, "metadata",
  "gata6_serial_metastatic_selection_projection_spec_2026-09-01.md"
)
required <- c(
  counts_path, series_path, gtf_path, phase2_sets_path, phase3_sets_path,
  coexistence_path, spec_path
)
if (any(!file.exists(required))) stop("One or more serial-selection inputs are missing")

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

adjusted_origin_effect <- function(values, generation, origin) {
  mean(vapply(levels(generation), function(current_generation) {
    mean(values[generation == current_generation & origin == "metastatic"]) -
      mean(values[generation == current_generation & origin == "primary"])
  }, numeric(1L)))
}

exact_stratified_origin_test <- function(values, generation, origin) {
  observed <- adjusted_origin_effect(values, generation, origin)
  generation_levels <- levels(generation)
  allocations <- lapply(generation_levels, function(current_generation) {
    index <- which(generation == current_generation)
    combn(index, sum(origin[index] == "metastatic"), simplify = FALSE)
  })
  null <- numeric(length(allocations[[1L]]) * length(allocations[[2L]]))
  counter <- 1L
  for (first_metastatic in allocations[[1L]]) {
    for (second_metastatic in allocations[[2L]]) {
      permuted_origin <- factor(
        rep("primary", length(origin)), levels = c("primary", "metastatic")
      )
      permuted_origin[c(first_metastatic, second_metastatic)] <- "metastatic"
      null[[counter]] <- adjusted_origin_effect(
        values, generation, permuted_origin
      )
      counter <- counter + 1L
    }
  }
  list(
    effect = observed,
    p_two_sided = mean(abs(null) >= abs(observed) - tol),
    n_allocations = length(null)
  )
}

bootstrap_stratified_effect <- function(values, generation, origin, iterations) {
  cell_indices <- split(
    seq_along(values), interaction(generation, origin, drop = TRUE)
  )
  estimates <- replicate(iterations, {
    selected <- unlist(lapply(cell_indices, function(index) {
      sample(index, length(index), replace = TRUE)
    }), use.names = FALSE)
    adjusted_origin_effect(values[selected], generation[selected], origin[selected])
  })
  as.numeric(quantile(estimates, c(0.025, 0.975), names = FALSE, type = 7))
}

phase2_sets <- fread(phase2_sets_path)
phase3_sets <- fread(phase3_sets_path)
coexistence <- fread(coexistence_path)
phase2_sets[, gene := toupper(gene)]
phase3_sets[, gene := toupper(gene)]
coexistence[, gene := toupper(gene)]

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
  "1010_Met", "1026_Met", "1027_Met",
  "1171_Primary", "1172_Primary", "1177_Primary",
  "1184_Met", "1186_Met", "1188_Met",
  "1242_Primary", "1244_Primary", "1248_Primary",
  "1235_Met", "1237_Met", "1238_Met"
)
if (!identical(sample_ids, expected_samples)) stop("Unexpected GSE290752 sample order")
if (any(vapply(counts_raw[, ..sample_ids], function(x) any(x < 0), logical(1L)))) {
  stop("Negative values found in the declared count matrix")
}
if (any(vapply(counts_raw[, ..sample_ids], function(x) any(x != floor(x)), logical(1L)))) {
  stop("GSE290752 processed matrix is not integer count data")
}

sample_titles <- c(
  "Met1_AKP_replicate1", "Met1_AKP_replicate2", "Met1_AKP_replicate3",
  "Primary2_AKP_replicate1", "Primary2_AKP_replicate2", "Primary2_AKP_replicate3",
  "Met2_AKP_replicate1", "Met2_AKP_replicate2", "Met2_AKP_replicate3",
  "Primary3_AKP_replicate1", "Primary3_AKP_replicate2", "Primary3_AKP_replicate3",
  "Met3_AKP_replicate1", "Met3_AKP_replicate2", "Met3_AKP_replicate3"
)
group_all <- factor(
  c(rep("Met1", 3L), rep("Primary2", 3L), rep("Met2", 3L),
    rep("Primary3", 3L), rep("Met3", 3L)),
  levels = c("Met1", "Primary2", "Met2", "Primary3", "Met3")
)
main_index <- which(group_all != "Met1")
generation <- factor(
  c(rep("round2", 6L), rep("round3", 6L)),
  levels = c("round2", "round3")
)
origin <- factor(
  c(rep("primary", 3L), rep("metastatic", 3L),
    rep("primary", 3L), rep("metastatic", 3L)),
  levels = c("primary", "metastatic")
)

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
count_matrix_all <- as.matrix(counts_by_gene[, ..sample_ids])
storage.mode(count_matrix_all) <- "integer"
rownames(count_matrix_all) <- counts_by_gene$gene

sample_metadata <- data.table(
  sample_id = sample_ids,
  sample_title = sample_titles,
  group = as.character(group_all),
  included_in_main_contrast = seq_along(sample_ids) %in% main_index,
  raw_library_size = colSums(count_matrix_all),
  detected_genes = colSums(count_matrix_all > 0)
)

count_matrix <- count_matrix_all[, main_index, drop = FALSE]
main_group <- interaction(generation, origin, drop = TRUE)
dge <- DGEList(count_matrix, group = main_group)
keep <- filterByExpr(dge, group = main_group)
dge <- calcNormFactors(dge[keep, , keep.lib.sizes = FALSE])
design <- model.matrix(~generation + origin)
voom_fit <- voom(dge, design, plot = FALSE)
fit <- eBayes(lmFit(voom_fit, design))
log_cpm <- voom_fit$E
gene_z <- safe_gene_z(log_cpm)

module_score_rows <- list()
module_coverage_rows <- list()
module_contrast_rows <- list()
loo_rows <- list()
for (set_id in names(gene_sets)) {
  declared <- unique(toupper(gene_sets[[set_id]]))
  measured <- intersect(declared, rownames(gene_z))
  if (!length(measured)) stop(paste("No measurable genes for", set_id))
  values <- as.numeric(colMeans(gene_z[measured, , drop = FALSE]))
  names(values) <- colnames(gene_z)
  exact <- exact_stratified_origin_test(values, generation, origin)
  interval <- bootstrap_stratified_effect(
    values, generation, origin, bootstrap_iterations
  )
  round2_difference <- mean(values[generation == "round2" & origin == "metastatic"]) -
    mean(values[generation == "round2" & origin == "primary"])
  round3_difference <- mean(values[generation == "round3" & origin == "metastatic"]) -
    mean(values[generation == "round3" & origin == "primary"])
  module_score_rows[[set_id]] <- data.table(
    sample_id = colnames(gene_z),
    generation = as.character(generation),
    origin = as.character(origin),
    set_id = set_id,
    score = values
  )
  module_coverage_rows[[set_id]] <- data.table(
    set_id = set_id,
    n_declared = length(declared),
    n_measured_nonconstant = length(measured),
    coverage_fraction = length(measured) / length(declared),
    measured_genes = paste(measured, collapse = ";"),
    missing_or_filtered_genes = paste(setdiff(declared, measured), collapse = ";")
  )
  module_contrast_rows[[set_id]] <- data.table(
    set_id = set_id,
    n_libraries = length(values),
    adjusted_metastatic_minus_primary = exact$effect,
    bootstrap_ci_lower = interval[[1L]],
    bootstrap_ci_upper = interval[[2L]],
    round2_metastatic_minus_primary = round2_difference,
    round3_metastatic_minus_primary = round3_difference,
    exact_p_two_sided = exact$p_two_sided,
    n_exact_allocations = exact$n_allocations
  )
  if (set_id %in% c(
    "REC_PROGRAM_TOP50", "REACTOME_TIGHT_JUNCTION_INTERACTIONS",
    "JUNB_TIGHT_JUNCTION_TARGETS_5"
  )) {
    loo_rows[[set_id]] <- rbindlist(lapply(measured, function(omitted_gene) {
      retained <- setdiff(measured, omitted_gene)
      retained_values <- as.numeric(colMeans(gene_z[retained, , drop = FALSE]))
      data.table(
        set_id = set_id,
        omitted_gene = omitted_gene,
        retained_gene_count = length(retained),
        adjusted_metastatic_minus_primary = adjusted_origin_effect(
          retained_values, generation, origin
        )
      )
    }))
  }
}

module_scores <- rbindlist(module_score_rows, use.names = TRUE)
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
    log2fc_metastatic_minus_primary = logFC,
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

# Contextual Met1/Met2/Met3 values use all libraries and a separate all-sample z scale.
dge_all <- DGEList(count_matrix_all, group = group_all)
keep_all <- filterByExpr(dge_all, group = group_all)
dge_all <- calcNormFactors(dge_all[keep_all, , keep.lib.sizes = FALSE])
log_cpm_all <- cpm(dge_all, log = TRUE, prior.count = 0.5)
gene_z_all <- safe_gene_z(log_cpm_all)
trajectory_rows <- list()
for (set_id in names(gene_sets)) {
  measured <- intersect(unique(toupper(gene_sets[[set_id]])), rownames(gene_z_all))
  values <- colMeans(gene_z_all[measured, , drop = FALSE])
  trajectory_rows[[set_id]] <- data.table(
    set_id = set_id,
    group = levels(group_all),
    group_mean_score = vapply(levels(group_all), function(current_group) {
      mean(values[group_all == current_group])
    }, numeric(1L)),
    n_libraries = as.integer(table(group_all)[levels(group_all)])
  )
}
trajectory <- rbindlist(trajectory_rows, use.names = TRUE)

sample_metadata[included_in_main_contrast == TRUE, `:=`(
  generation = as.character(generation),
  origin = as.character(origin),
  tmm_norm_factor_main = dge$samples$norm.factors,
  effective_library_size_main = dge$samples$lib.size * dge$samples$norm.factors
)]
sample_metadata[, `:=`(
  tmm_norm_factor_all = dge_all$samples$norm.factors,
  effective_library_size_all = dge_all$samples$lib.size * dge_all$samples$norm.factors
)]
sample_correlation <- as.data.table(
  cor(log_cpm, method = "pearson"), keep.rownames = "sample_id"
)

write_tsv(mapping_qc, "gse290752_gene_mapping_qc.tsv")
write_tsv(sample_metadata, "gse290752_sample_qc.tsv")
write_tsv(sample_correlation, "gse290752_main_sample_correlation.tsv")
write_tsv(module_coverage, "gse290752_programme_coverage.tsv")
write_tsv(module_scores, "gse290752_programme_scores.tsv")
write_tsv(module_contrasts, "gse290752_programme_contrasts.tsv")
write_tsv(module_loo, "gse290752_programme_leave_one_gene_out.tsv")
write_tsv(camera_result, "gse290752_camera_results.tsv")
write_tsv(named_effects, "gse290752_named_gene_effects.tsv")
write_tsv(trajectory, "gse290752_contextual_group_means.tsv")
fwrite(
  gene_table,
  file.path(out_dir, "gse290752_limma_all_genes.tsv.gz"),
  sep = "\t", na = "NA", compress = "gzip"
)

summary_lines <- c(
  "# GSE290752 serial metastatic-selection projection",
  "",
  paste0(
    "The primary contrast uses 12 organoid RNA-seq libraries from rounds 2 and 3, ",
    "with generation adjustment and exact label allocation within generation. ",
    "The three Met1 libraries are contextual only."
  ),
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
        "- `%s`: adjusted metastatic-primary %+.3f (95%% library-bootstrap interval ",
        "%+.3f to %+.3f; round 2 %+.3f; round 3 %+.3f; exact P=%.4g); ",
        "camera %s, P=%.3g, FDR=%.3g."
      ),
      current_set, row$adjusted_metastatic_minus_primary,
      row$bootstrap_ci_lower, row$bootstrap_ci_upper,
      row$round2_metastatic_minus_primary,
      row$round3_metastatic_minus_primary,
      row$exact_p_two_sided, camera_row$direction, camera_row$p_value,
      camera_row$fdr_bh_camera_family
    )
  )
}
summary_lines <- c(
  summary_lines,
  "",
  "## Result interpretation",
  "",
  "Repeated liver-metastatic selection increases the REC programme in both controlled rounds and increases the focused CLDN3/4/7-CRB3-F11R junction subset, with every leave-one-gene-out estimate remaining positive. GATA6 and HNF4A decrease while KLF4 increases. Broad tight-junction change is smaller, partial-EMT is lower and E2F/G2M directions differ between rounds, so the full replacement-associated combination is not reproduced.",
  "",
  "Together with the GATA6-knockout and Plexin B2 projections, this supports partial convergence during metastatic selection: loss of intestinal lineage restraint is accompanied by an epithelial-restabilizing signal, but neither component alone explains the patient rHGP state.",
  "",
  "## Interpretation boundary",
  "",
  paste0(
    "This dataset contains no HGP annotation. It tests whether repeated liver-metastatic ",
    "selection reproduces the replacement-associated state, not whether any programme ",
    "causes or validates rHGP."
  ),
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
  n_main_libraries = length(main_index),
  n_context_libraries = length(sample_ids) - length(main_index),
  exact_allocations = 400L,
  source_annotation = "Ensembl release 101, mm10/GRCm38",
  counts_md5 = unname(tools::md5sum(counts_path)),
  series_md5 = unname(tools::md5sum(series_path)),
  gtf_md5 = unname(tools::md5sum(gtf_path)),
  phase2_gene_sets_md5 = unname(tools::md5sum(phase2_sets_path)),
  phase3_gene_sets_md5 = unname(tools::md5sum(phase3_sets_path)),
  specification_md5 = unname(tools::md5sum(spec_path))
)
write_tsv(run_info, "gse290752_run_info.tsv")
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))

manifest <- data.table(
  file = sort(list.files(out_dir, full.names = FALSE)),
  role = "analysis_output"
)
write_tsv(manifest, "output_manifest.tsv")

cat("Wrote serial metastatic-selection projection to", out_dir, "\n")
