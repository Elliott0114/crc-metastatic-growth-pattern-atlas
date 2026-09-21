#!/usr/bin/env Rscript

# Reader-facing specificity and sensitivity analyses for the REC programme.
# The main comparisons are deliberately simple: state-program projections,
# an epithelial-content covariate, component scores, and scale-free ranks.

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(Matrix)
  library(readxl)
})

options(stringsAsFactors = FALSE)
set.seed(42)

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
ROOT <- if (length(file_arg)) {
  normalizePath(file.path(dirname(sub("^--file=", "", file_arg[[1]])), ".."), mustWork = TRUE)
} else {
  normalizePath(".", mustWork = TRUE)
}

OUT <- file.path(ROOT, "analysis_results", "rec_program_reader_facing_specificity")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(x, filename) {
  fwrite(as.data.table(x), file.path(OUT, filename), sep = "\t", na = "NA")
}

write_tsv_gz <- function(x, filename) {
  connection <- gzfile(
    file.path(OUT, filename), open = "wt", encoding = "UTF-8"
  )
  on.exit(close(connection), add = TRUE)
  write.table(
    as.data.frame(x), connection, sep = "\t", quote = FALSE,
    row.names = FALSE, na = "NA"
  )
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
  z <- sweep(z, 1L, sds[keep], "/")
  z
}

cosine_similarity <- function(x, y) {
  sum(x * y) / sqrt(sum(x^2) * sum(y^2))
}

rank_states <- function(iss_values, ogden_values, metric = "pearson") {
  scores <- vapply(seq_len(nrow(ogden_values)), function(i) {
    if (metric == "pearson") {
      cor(iss_values, ogden_values[i, ], method = "pearson")
    } else if (metric == "spearman") {
      cor(iss_values, ogden_values[i, ], method = "spearman")
    } else {
      cosine_similarity(iss_values, ogden_values[i, ])
    }
  }, numeric(1))
  data.table(
    Ogden_state = rownames(ogden_values),
    similarity = scores,
    rank = frank(-scores, ties.method = "min")
  )[order(rank, Ogden_state)]
}

exact_two_group <- function(values, labels, positive_label = "rHGP") {
  values <- as.numeric(values)
  labels <- as.character(labels)
  n_positive <- sum(labels == positive_label)
  observed <- mean(values[labels == positive_label]) - mean(values[labels != positive_label])
  allocations <- combn(seq_along(values), n_positive)
  positive_sums <- colSums(matrix(values[allocations], nrow = n_positive))
  null <- positive_sums / n_positive -
    (sum(values) - positive_sums) / (length(values) - n_positive)
  data.table(
    effect = observed,
    exact_p_one_sided = mean(null >= observed - sqrt(.Machine$double.eps)),
    exact_p_two_sided = mean(abs(null) >= abs(observed) - sqrt(.Machine$double.eps)),
    n_assignments = length(null)
  )
}

exact_paired <- function(differences) {
  differences <- as.numeric(differences)
  signs <- as.matrix(expand.grid(rep(list(c(-1, 1)), length(differences))))
  null <- rowMeans(sweep(signs, 2L, differences, "*"))
  observed <- mean(differences)
  data.table(
    effect = observed,
    exact_p_one_sided = mean(null >= observed - sqrt(.Machine$double.eps)),
    exact_p_two_sided = mean(abs(null) >= abs(observed) - sqrt(.Machine$double.eps)),
    n_assignments = length(null)
  )
}

exact_adjusted_unpaired <- function(outcome, covariate, labels, positive_label = "rHGP") {
  outcome <- as.numeric(outcome)
  covariate <- as.numeric(covariate)
  group <- as.integer(as.character(labels) == positive_label)
  observed <- unname(coef(lm(outcome ~ group + covariate))["group"])
  allocations <- combn(seq_along(group), sum(group == 1L))
  null <- apply(allocations, 2L, function(index) {
    permuted_group <- integer(length(group))
    permuted_group[index] <- 1L
    unname(coef(lm(outcome ~ permuted_group + covariate))["permuted_group"])
  })
  data.table(
    effect = observed,
    exact_p_two_sided = mean(abs(null) >= abs(observed) - sqrt(.Machine$double.eps)),
    n_assignments = length(null)
  )
}

exact_adjusted_paired <- function(outcome_difference, covariate_difference) {
  outcome_difference <- as.numeric(outcome_difference)
  covariate_difference <- as.numeric(covariate_difference)
  observed <- unname(coef(lm(outcome_difference ~ covariate_difference))["(Intercept)"])
  signs <- as.matrix(expand.grid(rep(list(c(-1, 1)), length(outcome_difference))))
  null <- apply(signs, 1L, function(sign_vector) {
    permuted_outcome <- outcome_difference * sign_vector
    permuted_covariate <- covariate_difference * sign_vector
    unname(coef(lm(permuted_outcome ~ permuted_covariate))["(Intercept)"])
  })
  data.table(
    effect = observed,
    exact_p_two_sided = mean(abs(null) >= abs(observed) - sqrt(.Machine$double.eps)),
    n_assignments = length(null)
  )
}

score_unpaired_program <- function(gene_z, genes, metadata, programme, family) {
  matched <- intersect(unique(toupper(genes)), rownames(gene_z))
  if (length(matched) < 3L) return(NULL)
  scores <- colMeans(gene_z[matched, , drop = FALSE], na.rm = TRUE)
  test <- exact_two_group(scores[metadata$sample], metadata$hgp, "rHGP")
  data.table(
    dataset = "GSE151165 bulk tumour",
    programme = programme,
    family = family,
    n_defined = length(unique(genes)),
    n_measured = length(matched),
    effect = test$effect,
    exact_p_one_sided = test$exact_p_one_sided,
    exact_p_two_sided = test$exact_p_two_sided,
    n_units = nrow(metadata)
  )
}

score_paired_program <- function(gene_z, genes, sample_metadata, programme, family) {
  matched <- intersect(unique(toupper(genes)), rownames(gene_z))
  if (length(matched) < 3L) return(NULL)
  scores <- colMeans(gene_z[matched, , drop = FALSE], na.rm = TRUE)
  score_dt <- copy(sample_metadata)
  score_dt[, score := as.numeric(scores[sample_id])]
  paired <- dcast(score_dt, patient ~ region, value.var = "score")
  paired <- paired[complete.cases(micro_tumour, macro_tumour)]
  test <- exact_paired(paired$macro_tumour - paired$micro_tumour)
  data.table(
    dataset = "GSE294385 paired lesion regions",
    programme = programme,
    family = family,
    n_defined = length(unique(genes)),
    n_measured = length(matched),
    effect = test$effect,
    exact_p_one_sided = test$exact_p_one_sided,
    exact_p_two_sided = test$exact_p_two_sided,
    n_units = nrow(paired)
  )
}

# -----------------------------------------------------------------------------
# 1. NC3-to-REC mapping without relying on one metric or the six source markers
# -----------------------------------------------------------------------------

centroids <- fread(
  cmd = paste("zcat", shQuote(file.path(
    ROOT, "analysis_results", "cross_platform_state_anchor", "state_centroids_long.tsv.gz"
  )))
)
shared_audit <- fread(file.path(
  ROOT, "analysis_results", "cross_platform_state_anchor", "shared_gene_audit.tsv"
))
selected_genes <- shared_audit[tolower(as.character(selected_for_primary_anchor)) == "true", gene]
centroids <- centroids[gene %in% selected_genes]

centroid_matrix <- function(dataset_name) {
  x <- dcast(
    centroids[dataset == dataset_name], state ~ gene,
    value.var = "patient_balanced_log1p_CPM"
  )
  states <- x$state
  mat <- as.matrix(x[, -"state"])
  mode(mat) <- "numeric"
  rownames(mat) <- states
  gene_means <- colMeans(mat)
  gene_sds <- sqrt(colMeans(sweep(mat, 2L, gene_means, "-")^2))
  keep <- is.finite(gene_sds) & gene_sds > 0
  z <- sweep(mat[, keep, drop = FALSE], 2L, gene_means[keep], "-")
  z <- sweep(z, 2L, gene_sds[keep], "/")
  z
}

iss_z <- centroid_matrix("ISS")
ogden_z <- centroid_matrix("Ogden")
common_anchor_genes <- intersect(colnames(iss_z), colnames(ogden_z))
iss_z <- iss_z[, common_anchor_genes, drop = FALSE]
ogden_z <- ogden_z[, common_anchor_genes, drop = FALSE]

metric_rows <- rbindlist(lapply(c("pearson", "spearman", "cosine"), function(metric) {
  ranked <- rank_states(iss_z["NC3", ], ogden_z, metric)
  rec <- ranked[Ogden_state == "REC"]
  second <- ranked[rank == 2L][1]
  data.table(
    analysis = metric,
    n_genes = ncol(iss_z),
    rec_similarity = rec$similarity,
    rec_rank = rec$rank,
    second_state = second$Ogden_state,
    second_similarity = second$similarity,
    rec_minus_second = rec$similarity - second$similarity
  )
}))

source_markers <- c("KRT18", "KRT19", "KRT20", "SLPI", "CDX2", "ASCL2")
marker_removed_genes <- setdiff(common_anchor_genes, source_markers)
removed_rank <- rank_states(
  iss_z["NC3", marker_removed_genes],
  ogden_z[, marker_removed_genes, drop = FALSE],
  "pearson"
)
removed_rec <- removed_rank[Ogden_state == "REC"]
removed_second <- removed_rank[rank == 2L][1]
metric_rows <- rbind(
  metric_rows,
  data.table(
    analysis = "pearson_without_six_source_markers",
    n_genes = length(marker_removed_genes),
    rec_similarity = removed_rec$similarity,
    rec_rank = removed_rec$rank,
    second_state = removed_second$Ogden_state,
    second_similarity = removed_second$similarity,
    rec_minus_second = removed_rec$similarity - removed_second$similarity
  )
)
write_tsv(metric_rows, "state_mapping_metric_and_marker_sensitivity.tsv")

loo_gene <- rbindlist(lapply(common_anchor_genes, function(gene) {
  keep <- setdiff(common_anchor_genes, gene)
  ranked <- rank_states(iss_z["NC3", keep], ogden_z[, keep, drop = FALSE], "pearson")
  data.table(
    omitted_gene = gene,
    top_state = ranked$Ogden_state[[1]],
    top_similarity = ranked$similarity[[1]],
    second_state = ranked$Ogden_state[[2]],
    top_minus_second = ranked$similarity[[1]] - ranked$similarity[[2]]
  )
}))
write_tsv(loo_gene, "state_mapping_leave_one_gene_out.tsv")

bootstrap_n <- 10000L
bootstrap_rows <- rbindlist(lapply(seq_len(bootstrap_n), function(iteration) {
  sampled <- sample(common_anchor_genes, length(common_anchor_genes), replace = TRUE)
  ranked <- rank_states(iss_z["NC3", sampled], ogden_z[, sampled, drop = FALSE], "pearson")
  data.table(
    iteration = iteration,
    top_state = ranked$Ogden_state[[1]],
    rec_similarity = ranked[Ogden_state == "REC", similarity],
    top_minus_second = ranked$similarity[[1]] - ranked$similarity[[2]]
  )
}))
write_tsv_gz(bootstrap_rows, "state_mapping_gene_bootstrap.tsv.gz")
bootstrap_summary <- data.table(
  n_iterations = bootstrap_n,
  rec_top_fraction = mean(bootstrap_rows$top_state == "REC"),
  rec_similarity_median = median(bootstrap_rows$rec_similarity),
  rec_similarity_q025 = quantile(bootstrap_rows$rec_similarity, 0.025),
  rec_similarity_q975 = quantile(bootstrap_rows$rec_similarity, 0.975),
  top_margin_median = median(bootstrap_rows$top_minus_second),
  positive_margin_fraction = mean(bootstrap_rows$top_minus_second > 0)
)
write_tsv(bootstrap_summary, "state_mapping_gene_bootstrap_summary.tsv")

# -----------------------------------------------------------------------------
# 2. Prepare GSE151165 tumour expression
# -----------------------------------------------------------------------------

gse151_path <- file.path(
  ROOT, "data_sources", "GSE151165", "GSE151165_RNA_seq.raw_read_count.xlsx"
)
raw <- as.data.table(read_excel(gse151_path, .name_repair = "unique_quiet"))
sample_columns <- grep("^KR-", names(raw), value = TRUE)
setnames(raw, 1L, "source_gene")
raw[, gene := toupper(trimws(as.character(source_gene)))]
raw <- raw[!is.na(gene) & gene != ""]
bulk_counts_dt <- raw[, lapply(.SD, function(x) sum(as.numeric(x), na.rm = TRUE)),
                      by = gene, .SDcols = sample_columns]
bulk_counts <- as.matrix(bulk_counts_dt[, ..sample_columns])
mode(bulk_counts) <- "numeric"
rownames(bulk_counts) <- bulk_counts_dt$gene

bulk_metadata_all <- data.table(
  sample = sample_columns,
  class = c(rep("D_N", 9L), rep("D_T", 9L), rep("R_N", 6L), rep("R_T", 6L))
)
bulk_metadata_all[, hgp := fifelse(startsWith(class, "R"), "rHGP", "dHGP")]
bulk_metadata_all[, tissue := fifelse(endsWith(class, "_T"), "tumour", "adjacent_liver")]
bulk_metadata <- bulk_metadata_all[tissue == "tumour"]
bulk_counts_tumour <- bulk_counts[, bulk_metadata$sample, drop = FALSE]

y_bulk <- DGEList(bulk_counts_tumour, group = bulk_metadata$hgp)
keep_bulk <- filterByExpr(y_bulk, group = bulk_metadata$hgp)
y_bulk <- calcNormFactors(y_bulk[keep_bulk, , keep.lib.sizes = FALSE])
bulk_effective_library <- y_bulk$samples$lib.size * y_bulk$samples$norm.factors
names(bulk_effective_library) <- rownames(y_bulk$samples)
bulk_log_cpm <- cpm(
  bulk_counts_tumour,
  lib.size = bulk_effective_library[colnames(bulk_counts_tumour)],
  log = TRUE,
  prior.count = 1
)
bulk_gene_z <- safe_z_rows(bulk_log_cpm, population_sd = FALSE)

# -----------------------------------------------------------------------------
# 3. Aggregate all genes in GSE294385 once, then reuse the patient-level matrix
# -----------------------------------------------------------------------------

spatial_cache <- file.path(OUT, "gse294385_patient_region_all_gene_counts.tsv.gz")
spatial_meta_cache <- file.path(OUT, "gse294385_patient_region_metadata.tsv")

if (!file.exists(spatial_cache) || !file.exists(spatial_meta_cache)) {
  message("Aggregating all GSE294385 genes to patient-region pseudobulks")
  manifest <- fread(file.path(ROOT, "metadata", "gse294385_liver_sample_manifest.tsv"))
  manifest <- manifest[tolower(as.character(selected_for_paired_extension)) == "true"]
  annotation_path <- file.path(
    ROOT, "data_sources", "Liu_2026_GSE294385", "visium_liver_meta_after_qc.tsv.gz"
  )
  annotation <- fread(cmd = paste("zcat", shQuote(annotation_path)))
  annotation <- annotation[Layer3 %in% c(
    "Liver micrometastasis tumor", "Liver macrometastasis tumor"
  )]
  annotation[, region := fifelse(
    Layer3 == "Liver micrometastasis tumor", "micro_tumour", "macro_tumour"
  )]

  sample_rows <- vector("list", nrow(manifest) * 2L)
  row_index <- 0L
  for (i in seq_len(nrow(manifest))) {
    sample_id <- manifest$sample[[i]]
    patient_id <- manifest$patient[[i]]
    matrix_dir <- file.path(
      ROOT, "data_sources", "Liu_2026_GSE294385", "extracted",
      sample_id, "filtered_feature_bc_matrix"
    )
    barcodes <- scan(
      gzfile(file.path(matrix_dir, "barcodes.tsv.gz")), what = character(),
      quiet = TRUE, encoding = "UTF-8"
    )
    feature_path <- file.path(matrix_dir, "features.tsv.gz")
    features <- fread(cmd = paste("zcat", shQuote(feature_path)), header = FALSE)
    genes <- toupper(as.character(features[[2]]))
    mat <- readMM(gzfile(file.path(matrix_dir, "matrix.mtx.gz")))
    if (nrow(mat) != length(genes) || ncol(mat) != length(barcodes)) {
      stop("GSE294385 matrix dimensions do not match features/barcodes for ", sample_id)
    }
    sample_annotation <- annotation[sample == sample_id]
    region_by_barcode <- sample_annotation$region[match(barcodes, sample_annotation$spot_barcode)]
    for (region_name in c("micro_tumour", "macro_tumour")) {
      selected_columns <- which(region_by_barcode == region_name)
      if (!length(selected_columns)) next
      region_counts <- Matrix::rowSums(mat[, selected_columns, drop = FALSE])
      dt <- data.table(gene = genes, raw_count = as.numeric(region_counts))[
        !is.na(gene) & gene != "", .(raw_count = sum(raw_count)), by = gene
      ]
      dt[, `:=`(
        sample = sample_id,
        patient = patient_id,
        region = region_name,
        n_spots = length(selected_columns),
        library_size = sum(raw_count)
      )]
      row_index <- row_index + 1L
      sample_rows[[row_index]] <- dt
    }
    rm(mat)
    invisible(gc())
  }
  sample_long <- rbindlist(sample_rows[seq_len(row_index)], use.names = TRUE)
  sample_meta <- unique(sample_long[, .(sample, patient, region, n_spots, library_size)])
  patient_long <- sample_long[, .(raw_count = sum(raw_count)), by = .(patient, region, gene)]
  patient_meta <- sample_meta[, .(
    n_spots = sum(n_spots),
    library_size = sum(library_size)
  ), by = .(patient, region)]
  eligible_patients <- patient_meta[, .(n_regions = uniqueN(region)), by = patient][n_regions == 2L, patient]
  patient_long <- patient_long[patient %in% eligible_patients]
  patient_meta <- patient_meta[patient %in% eligible_patients]
  write_tsv_gz(patient_long, basename(spatial_cache))
  write_tsv(patient_meta, basename(spatial_meta_cache))
} else {
  patient_long <- fread(cmd = paste("zcat", shQuote(spatial_cache)))
  patient_meta <- fread(spatial_meta_cache)
}

patient_long[, sample_id := paste(patient, region, sep = "__")]
patient_meta[, sample_id := paste(patient, region, sep = "__")]
spatial_wide <- dcast(patient_long, gene ~ sample_id, value.var = "raw_count", fill = 0)
spatial_genes <- spatial_wide$gene
spatial_counts <- as.matrix(spatial_wide[, -"gene"])
mode(spatial_counts) <- "numeric"
rownames(spatial_counts) <- spatial_genes
spatial_meta <- unique(patient_meta[, .(patient, region, sample_id, library_size, n_spots)])
spatial_meta <- spatial_meta[match(colnames(spatial_counts), sample_id)]

# The original locked projection required a gene to be present in every source
# feature table. Apply the same platform-coverage rule to all new programmes.
spatial_manifest <- fread(file.path(
  ROOT, "metadata", "gse294385_liver_sample_manifest.tsv"
))
spatial_manifest <- spatial_manifest[
  tolower(as.character(selected_for_paired_extension)) == "true"
]
spatial_feature_sets <- lapply(spatial_manifest$sample, function(sample_id) {
  feature_path <- file.path(
    ROOT, "data_sources", "Liu_2026_GSE294385", "extracted",
    sample_id, "filtered_feature_bc_matrix", "features.tsv.gz"
  )
  unique(toupper(fread(
    cmd = paste("zcat", shQuote(feature_path)), header = FALSE,
    select = 2L
  )[[1]]))
})
common_spatial_features <- Reduce(intersect, spatial_feature_sets)
spatial_counts <- spatial_counts[
  intersect(rownames(spatial_counts), common_spatial_features), , drop = FALSE
]
spatial_log_cpm <- log2(sweep(spatial_counts + 0.5, 2L, spatial_meta$library_size + 1, "/") * 1e6)
spatial_gene_z <- safe_z_rows(spatial_log_cpm, population_sd = TRUE)

# Confirm that the new full-gene aggregation reproduces the locked REC result.
rec_program <- fread(file.path(
  ROOT, "analysis_results", "ogden_anchor_program", "selected_anchor_program_top50.tsv"
))[order(program_rank)]
rec_program[, gene := toupper(gene)]
locked_spatial <- fread(file.path(
  ROOT, "analysis_results", "gse294385_rec_program_extension", "paired_contrasts.tsv"
))[comparison == "macro_minus_micro_tumour" & program == "rec_equal_gene"]
recomputed_spatial <- score_paired_program(
  spatial_gene_z, rec_program$gene, spatial_meta, "Fixed REC 50", "primary"
)
if (abs(recomputed_spatial$effect - locked_spatial$mean_difference) > 1e-8) {
  stop(sprintf(
    "Full-gene spatial aggregation does not reproduce the locked REC effect: %.12f vs %.12f",
    recomputed_spatial$effect, locked_spatial$mean_difference
  ))
}

# Programme length is an operational choice. Test whether the two principal
# external directions depend on choosing exactly 50 genes, without changing the
# frozen 50-gene endpoint used elsewhere in the article.
eligible_programme <- fread(file.path(
  ROOT, "analysis_results", "ogden_anchor_program", "eligible_program_candidates.tsv"
))[order(-t, -logFC)]
eligible_programme[, gene := toupper(gene)]
if (!identical(eligible_programme$gene[seq_len(50L)], rec_program$gene)) {
  stop("The first 50 eligible REC genes do not reproduce the frozen programme")
}
programme_lengths <- c(25L, 50L, 75L, 100L)
programme_length_effects <- rbindlist(lapply(programme_lengths, function(n_genes) {
  genes <- eligible_programme$gene[seq_len(n_genes)]
  rbind(
    score_unpaired_program(
      bulk_gene_z, genes, bulk_metadata,
      sprintf("Top %d REC genes", n_genes), "programme length"
    ),
    score_paired_program(
      spatial_gene_z, genes, spatial_meta,
      sprintf("Top %d REC genes", n_genes), "programme length"
    )
  )[, programme_size := n_genes]
}))
setcolorder(
  programme_length_effects,
  c("dataset", "programme", "programme_size", "family", "n_defined",
    "n_measured", "effect", "exact_p_one_sided", "exact_p_two_sided", "n_units")
)
write_tsv(programme_length_effects, "rec_program_length_sensitivity.tsv")
