#!/usr/bin/env Rscript

# Place the fixed REC programme in tumour-versus-adjacent-liver context in GSE151165.

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
out_dir <- file.path(root, "analysis_results", "rec_program_gse151165_tissue_context")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

bootstrap_replicates <- 10000L

write_tsv <- function(x, filename) {
  fwrite(as.data.table(x), file.path(out_dir, filename), sep = "\t", na = "NA")
}

exact_two_group <- function(values, labels, positive_label) {
  values <- as.numeric(values)
  labels <- as.character(labels)
  n_positive <- sum(labels == positive_label)
  observed <- mean(values[labels == positive_label]) - mean(values[labels != positive_label])
  allocations <- utils::combn(seq_along(values), n_positive)
  positive_sums <- colSums(matrix(values[allocations], nrow = n_positive))
  null_effects <- positive_sums / n_positive -
    (sum(values) - positive_sums) / (length(values) - n_positive)
  data.table(
    effect = observed,
    exact_p_two_sided = mean(abs(null_effects) >= abs(observed) - sqrt(.Machine$double.eps)),
    n_allocations = length(null_effects)
  )
}

bootstrap_four_group <- function(dt, statistic, replicates) {
  group_rows <- split(seq_len(nrow(dt)), interaction(dt$tissue, dt$hgp, drop = TRUE))
  values <- replicate(replicates, {
    sampled_rows <- unlist(lapply(group_rows, function(rows) sample(rows, length(rows), replace = TRUE)))
    statistic(dt[sampled_rows])
  })
  as.numeric(quantile(values, c(0.025, 0.975), names = FALSE, type = 7))
}

raw <- as.data.table(read_excel(raw_path, .name_repair = "unique_quiet"))
sample_columns <- grep("^KR-", names(raw), value = TRUE)
if (length(sample_columns) != 30L) {
  stop("Expected 30 GSE151165 sample columns")
}

setnames(raw, 1L, "source_gene")
raw[, gene := toupper(trimws(as.character(source_gene)))]
raw <- raw[!is.na(gene) & gene != ""]
counts_table <- raw[
  , lapply(.SD, function(x) sum(as.numeric(x), na.rm = TRUE)),
  by = gene,
  .SDcols = sample_columns
]
counts <- as.matrix(counts_table[, ..sample_columns])
mode(counts) <- "numeric"
rownames(counts) <- counts_table$gene

metadata <- data.table(
  sample = sample_columns,
  class = c(rep("D_N", 9L), rep("D_T", 9L), rep("R_N", 6L), rep("R_T", 6L))
)
metadata[, hgp := fifelse(startsWith(class, "R"), "rHGP", "dHGP")]
metadata[, tissue := fifelse(endsWith(class, "_T"), "Tumour", "Adjacent liver")]

y <- DGEList(counts)
y <- calcNormFactors(y)
effective_library <- y$samples$lib.size * y$samples$norm.factors
names(effective_library) <- rownames(y$samples)
log_cpm <- cpm(counts, lib.size = effective_library[colnames(counts)], log = TRUE, prior.count = 1)

program <- fread(program_path)
program[, gene := toupper(gene)]
program_genes <- intersect(program$gene, rownames(log_cpm))
if (length(program_genes) < 45L) {
  stop("REC programme coverage unexpectedly below 45 genes")
}

gene_z <- t(scale(t(log_cpm[program_genes, , drop = FALSE])))
metadata[, rec_program_score := colMeans(gene_z, na.rm = TRUE)]
pooled_counts <- colSums(counts[program_genes, , drop = FALSE])
metadata[, pooled_program_log2_cpm := as.numeric(log2(
  (pooled_counts[sample] + 1) / effective_library[sample] * 1e6
))]

group_summary <- metadata[
  , .(
    n_samples = .N,
    mean_rec_score = mean(rec_program_score),
    sd_rec_score = sd(rec_program_score),
    median_rec_score = median(rec_program_score),
    mean_pooled_log2_cpm = mean(pooled_program_log2_cpm)
  ),
  by = .(tissue, hgp)
]

score_interaction <- function(dt) {
  means <- dcast(dt, tissue ~ hgp, value.var = "rec_program_score", fun.aggregate = mean)
  tumour_effect <- means[tissue == "Tumour", rHGP - dHGP]
  liver_effect <- means[tissue == "Adjacent liver", rHGP - dHGP]
  tumour_effect - liver_effect
}

tumour_exact <- exact_two_group(
  metadata[tissue == "Tumour", rec_program_score],
  metadata[tissue == "Tumour", hgp],
  "rHGP"
)
liver_exact <- exact_two_group(
  metadata[tissue == "Adjacent liver", rec_program_score],
  metadata[tissue == "Adjacent liver", hgp],
  "rHGP"
)

fit <- lm(rec_program_score ~ tissue * hgp, data = metadata)
fit_table <- as.data.table(coef(summary(fit)), keep.rownames = "term")
setnames(fit_table, c("Estimate", "Std. Error", "t value", "Pr(>|t|)"),
         c("estimate", "standard_error", "t_value", "p_value"))
interaction_row <- fit_table[grepl(":", term)]
interaction_ci <- bootstrap_four_group(metadata, score_interaction, bootstrap_replicates)

contrast_summary <- rbindlist(list(
  data.table(
    contrast = "rHGP_minus_dHGP_in_tumour",
    estimate = tumour_exact$effect,
    ci_lower = NA_real_,
    ci_upper = NA_real_,
    p_value = tumour_exact$exact_p_two_sided,
    test = "exact label permutation",
    n_assignments = tumour_exact$n_allocations
  ),
  data.table(
    contrast = "rHGP_minus_dHGP_in_adjacent_liver",
    estimate = liver_exact$effect,
    ci_lower = NA_real_,
    ci_upper = NA_real_,
    p_value = liver_exact$exact_p_two_sided,
    test = "exact label permutation",
    n_assignments = liver_exact$n_allocations
  ),
  data.table(
    contrast = "tumour_specific_HGP_amplification",
    estimate = score_interaction(metadata),
    ci_lower = interaction_ci[[1]],
    ci_upper = interaction_ci[[2]],
    p_value = interaction_row$p_value,
    test = "ordinary least-squares tissue-by-HGP interaction",
    n_assignments = NA_integer_
  )
), use.names = TRUE, fill = TRUE)

gene_effects <- rbindlist(lapply(program_genes, function(gene_symbol) {
  dt <- copy(metadata)
  dt[, value := as.numeric(log_cpm[gene_symbol, sample])]
  tumour_effect <- dt[tissue == "Tumour", mean(value[hgp == "rHGP"]) - mean(value[hgp == "dHGP"])]
  liver_effect <- dt[tissue == "Adjacent liver", mean(value[hgp == "rHGP"]) - mean(value[hgp == "dHGP"])]
  data.table(
    gene = gene_symbol,
    program_rank = program[gene == gene_symbol, program_rank],
    tumour_rHGP_minus_dHGP = tumour_effect,
    adjacent_liver_rHGP_minus_dHGP = liver_effect,
    tumour_specific_amplification = tumour_effect - liver_effect
  )
}))

coverage <- data.table(
  n_defined = nrow(program),
  n_matched = length(program_genes),
  coverage_fraction = length(program_genes) / nrow(program),
  matched_genes = paste(program_genes, collapse = ","),
  missing_genes = paste(setdiff(program$gene, program_genes), collapse = ",")
)

write_tsv(coverage, "program_coverage.tsv")
write_tsv(metadata[order(tissue, hgp, sample)], "sample_scores.tsv")
write_tsv(group_summary[order(tissue, hgp)], "group_summary.tsv")
write_tsv(contrast_summary, "contrast_summary.tsv")
write_tsv(fit_table, "interaction_model.tsv")
write_tsv(gene_effects[order(program_rank)], "gene_tissue_effects.tsv")

run_info <- data.table(
  analysis_date = "2026-08-23",
  statistical_unit = "public sample",
  n_samples = nrow(metadata),
  n_tumour = metadata[tissue == "Tumour", .N],
  n_adjacent_liver = metadata[tissue == "Adjacent liver", .N],
  n_program_genes = length(program_genes),
  bootstrap_replicates = bootstrap_replicates,
  seed = 42L,
  r_version = R.version.string,
  edgeR_version = as.character(packageVersion("edgeR")),
  raw_file_md5 = unname(tools::md5sum(raw_path)),
  programme_file_md5 = unname(tools::md5sum(program_path))
)
write_tsv(run_info, "run_info.tsv")
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))

summary_lines <- c(
  "# REC programme in tumour and adjacent liver in GSE151165",
  "",
  sprintf(
    "The fixed programme retained %d/%d genes across 30 samples (9 dHGP and 6 rHGP samples per tissue compartment).",
    length(program_genes), nrow(program)
  ),
  "",
  sprintf(
    "The rHGP-minus-dHGP programme difference was %+.3f in tumour (exact P=%.4f) and %+.3f in adjacent liver (exact P=%.4f).",
    tumour_exact$effect, tumour_exact$exact_p_two_sided,
    liver_exact$effect, liver_exact$exact_p_two_sided
  ),
  sprintf(
    "The tumour-specific amplification was %+.3f (bootstrap 95%% CI %+.3f to %+.3f; interaction P=%.4f).",
    score_interaction(metadata), interaction_ci[[1]], interaction_ci[[2]], interaction_row$p_value
  ),
  "",
  "The samples are independent public columns rather than paired tumour-adjacent specimens. The interaction therefore localises the stronger HGP contrast to tumour samples but is not a within-patient tissue comparison."
)
writeLines(summary_lines, file.path(out_dir, "analysis_summary.md"), useBytes = TRUE)

manifest <- data.table(
  file = c(
    "program_coverage.tsv", "sample_scores.tsv", "group_summary.tsv",
    "contrast_summary.tsv", "interaction_model.tsv", "gene_tissue_effects.tsv",
    "run_info.tsv", "sessionInfo.txt", "analysis_summary.md"
  ),
  role = c(
    "fixed programme coverage", "sample-level scores", "four-group descriptive summary",
    "tissue-specific HGP contrasts", "tissue-by-HGP model", "gene-level tissue contrasts",
    "environment and source fingerprints", "R environment", "human-readable result summary"
  )
)
fwrite(manifest, file.path(out_dir, "_analysis_outputs.md"), sep = "\t")

