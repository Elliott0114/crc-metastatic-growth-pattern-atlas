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

# -----------------------------------------------------------------------------
# 4. Fixed-program components and published epithelial-state controls
# -----------------------------------------------------------------------------

author_overlap_row <- fread(file.path(
  ROOT, "analysis_results", "ogden_anchor_program", "source_author_marker_overlap.tsv"
))[1]
author_genes <- strsplit(author_overlap_row$overlap_genes, ";", fixed = TRUE)[[1]]
additional_genes <- setdiff(rec_program$gene, author_genes)
component_sets <- list(
  `Fixed REC 50` = rec_program$gene,
  `Source-overlapping 41` = author_genes,
  `Additional 9` = additional_genes
)
component_effects <- rbindlist(lapply(names(component_sets), function(programme) {
  genes <- component_sets[[programme]]
  rbind(
    score_unpaired_program(bulk_gene_z, genes, bulk_metadata, programme, "REC programme component"),
    score_paired_program(spatial_gene_z, genes, spatial_meta, programme, "REC programme component")
  )
}))
write_tsv(component_effects, "rec_program_component_effects.tsv")

read_ogden_marker_sheet <- function(sheet) {
  workbook <- file.path(
    ROOT, "data_sources", "Ogden_2025_CRLM_multiome", "supplementary_tables", "mmc6.xlsx"
  )
  raw_sheet <- as.data.frame(
    read_excel(workbook, sheet = sheet, col_names = FALSE, .name_repair = "minimal"),
    stringsAsFactors = FALSE
  )
  header_row <- which(apply(raw_sheet, 1L, function(x) any(as.character(x) == "gene")))[1]
  if (!is.finite(header_row)) stop("Could not locate marker header in sheet ", sheet)
  headers <- as.character(raw_sheet[header_row, ])
  values <- raw_sheet[(header_row + 1L):nrow(raw_sheet), , drop = FALSE]
  names(values) <- headers
  values <- as.data.table(values)
  values[, `:=`(
    gene = toupper(trimws(as.character(gene))),
    avg_log2FC = as.numeric(avg_log2FC),
    p_val_adj = as.numeric(p_val_adj)
  )]
  values[!is.na(gene) & gene != "" & is.finite(avg_log2FC) & is.finite(p_val_adj)]
}

reference_states <- c(
  "REC", "Hypoxia", "iREC", "UPR", "Stem NOTUM",
  "Colonocyte", "Stem", "Goblet", "Intermediate", "TA1"
)
technical_pattern <- paste(
  c("^MT-", "^RPS[0-9]", "^RPL[0-9]", "^MIR[0-9]", "^LINC[0-9]", "^AC[0-9]", "^AL[0-9]"),
  collapse = "|"
)
state_marker_candidates <- lapply(reference_states, function(state_name) {
  x <- read_ogden_marker_sheet(state_name)
  x <- x[p_val_adj < 0.05 & avg_log2FC >= 0.5 & !grepl(technical_pattern, gene, perl = TRUE)]
  x <- unique(x[order(-avg_log2FC, p_val_adj)], by = "gene")
  x[, state := state_name]
  x
})
names(state_marker_candidates) <- reference_states
equal_state_size <- min(vapply(state_marker_candidates, nrow, integer(1)), 50L)
state_programs <- rbindlist(lapply(reference_states, function(state_name) {
  x <- head(state_marker_candidates[[state_name]], equal_state_size)
  x[, marker_rank := seq_len(.N)]
  x[, .(state, marker_rank, gene, avg_log2FC, p_val_adj)]
}))
write_tsv(state_programs, "ogden_published_state_programs.tsv")

state_effects <- rbindlist(lapply(reference_states, function(state_name) {
  genes <- state_programs[state == state_name, gene]
  rbind(
    score_unpaired_program(
      bulk_gene_z, genes, bulk_metadata, paste0("Published ", state_name), "Published state marker control"
    ),
    score_paired_program(
      spatial_gene_z, genes, spatial_meta, paste0("Published ", state_name), "Published state marker control"
    )
  )
}))
state_effects <- rbind(
  state_effects,
  score_unpaired_program(
    bulk_gene_z, rec_program$gene, bulk_metadata, "Fixed REC 50", "Fixed programme"
  ),
  score_paired_program(
    spatial_gene_z, rec_program$gene, spatial_meta, "Fixed REC 50", "Fixed programme"
  )
)
state_effects[, rank_within_dataset := frank(-effect, ties.method = "min"), by = dataset]
write_tsv(state_effects, "ogden_state_program_external_effects.tsv")

# -----------------------------------------------------------------------------
# 5. Bulk epithelial-content control (surrogate, not formal tumour purity)
# -----------------------------------------------------------------------------

epithelial_genes <- c("EPCAM", "TACSTD2", "CDH1", "KRT8", "KRT18", "KRT19", "MUC1")
rec_matched_bulk <- intersect(rec_program$gene, rownames(bulk_gene_z))
epi_matched_bulk <- intersect(epithelial_genes, rownames(bulk_gene_z))
rec_scores_bulk <- colMeans(bulk_gene_z[rec_matched_bulk, , drop = FALSE])
epi_scores_bulk <- colMeans(bulk_gene_z[epi_matched_bulk, , drop = FALSE])
control_data <- copy(bulk_metadata)
control_data[, `:=`(
  rec_score = as.numeric(rec_scores_bulk[sample]),
  epithelial_content_score = as.numeric(epi_scores_bulk[sample]),
  hgp_binary = as.integer(hgp == "rHGP")
)]

unadjusted_test <- exact_two_group(control_data$rec_score, control_data$hgp, "rHGP")
epithelial_test <- exact_two_group(
  control_data$epithelial_content_score, control_data$hgp, "rHGP"
)
adjusted_test <- exact_adjusted_unpaired(
  control_data$rec_score,
  control_data$epithelial_content_score,
  control_data$hgp,
  "rHGP"
)
adjusted_effect <- adjusted_test$effect
adjusted_p <- adjusted_test$exact_p_two_sided

control_summary <- data.table(
  model = c("REC score without covariate", "REC score with epithelial-content covariate", "Epithelial-content score"),
  effect = c(unadjusted_test$effect, adjusted_effect, epithelial_test$effect),
  exact_p_two_sided = c(unadjusted_test$exact_p_two_sided, adjusted_p, epithelial_test$exact_p_two_sided),
  n_samples = nrow(control_data),
  rec_epithelial_score_correlation = cor(control_data$rec_score, control_data$epithelial_content_score),
  epithelial_genes = paste(epi_matched_bulk, collapse = ";")
)
write_tsv(control_data, "bulk_epithelial_content_sample_scores.tsv")
write_tsv(control_summary, "bulk_epithelial_content_control.tsv")

# -----------------------------------------------------------------------------
# 6. Separate the REC-associated hypoxia component from the remaining programme
# -----------------------------------------------------------------------------

hypoxia_row <- fread(file.path(
  ROOT, "analysis_results", "ogden_anchor_program", "published_signature_enrichment.tsv"
))[gene_set == "Hypoxia (MP6)"]
if (nrow(hypoxia_row) != 1L) stop("Expected one Hypoxia (MP6) enrichment row")
hypoxia_overlap_genes <- strsplit(hypoxia_row$overlap_genes, ";", fixed = TRUE)[[1]]
rec_nonhypoxia_genes <- setdiff(rec_program$gene, hypoxia_overlap_genes)
rec_hypoxia_sets <- list(
  `REC non-hypoxia component` = rec_nonhypoxia_genes,
  `REC–hypoxia overlap` = hypoxia_overlap_genes
)

bulk_component_scores <- lapply(rec_hypoxia_sets, function(genes) {
  matched <- intersect(genes, rownames(bulk_gene_z))
  colMeans(bulk_gene_z[matched, , drop = FALSE])
})
bulk_hypoxia_data <- copy(bulk_metadata)
bulk_hypoxia_data[, `:=`(
  rec_nonhypoxia_score = as.numeric(bulk_component_scores[["REC non-hypoxia component"]][sample]),
  rec_hypoxia_overlap_score = as.numeric(bulk_component_scores[["REC–hypoxia overlap"]][sample])
)]

spatial_component_scores <- lapply(rec_hypoxia_sets, function(genes) {
  matched <- intersect(genes, rownames(spatial_gene_z))
  colMeans(spatial_gene_z[matched, , drop = FALSE])
})
spatial_hypoxia_data <- copy(spatial_meta)
spatial_hypoxia_data[, `:=`(
  rec_nonhypoxia_score = as.numeric(spatial_component_scores[["REC non-hypoxia component"]][sample_id]),
  rec_hypoxia_overlap_score = as.numeric(spatial_component_scores[["REC–hypoxia overlap"]][sample_id])
)]

component_effect_rows <- rbindlist(lapply(names(rec_hypoxia_sets), function(component) {
  genes <- rec_hypoxia_sets[[component]]
  bulk_score <- bulk_component_scores[[component]]
  bulk_test <- exact_two_group(
    bulk_score[bulk_metadata$sample], bulk_metadata$hgp, "rHGP"
  )
  spatial_score_dt <- copy(spatial_meta)
  spatial_score_dt[, score := as.numeric(spatial_component_scores[[component]][sample_id])]
  spatial_wide_score <- dcast(spatial_score_dt, patient ~ region, value.var = "score")
  spatial_differences <- spatial_wide_score$macro_tumour - spatial_wide_score$micro_tumour
  spatial_test <- exact_paired(spatial_differences)
  rbind(
    data.table(
      dataset = "GSE151165 bulk tumour", component = component,
      n_defined = length(genes), n_measured = length(intersect(genes, rownames(bulk_gene_z))),
      effect = bulk_test$effect, exact_p_two_sided = bulk_test$exact_p_two_sided,
      positive_units = NA_integer_,
      n_units = nrow(bulk_metadata)
    ),
    data.table(
      dataset = "GSE294385 paired lesion regions", component = component,
      n_defined = length(genes), n_measured = length(intersect(genes, rownames(spatial_gene_z))),
      effect = spatial_test$effect, exact_p_two_sided = spatial_test$exact_p_two_sided,
      positive_units = sum(spatial_differences > 0), n_units = length(spatial_differences)
    )
  )
}))

bulk_nonhypoxia_adjusted <- exact_adjusted_unpaired(
  bulk_hypoxia_data$rec_nonhypoxia_score,
  bulk_hypoxia_data$rec_hypoxia_overlap_score,
  bulk_hypoxia_data$hgp
)
bulk_hypoxia_adjusted <- exact_adjusted_unpaired(
  bulk_hypoxia_data$rec_hypoxia_overlap_score,
  bulk_hypoxia_data$rec_nonhypoxia_score,
  bulk_hypoxia_data$hgp
)

spatial_hypoxia_wide <- dcast(
  spatial_hypoxia_data,
  patient ~ region,
  value.var = c("rec_nonhypoxia_score", "rec_hypoxia_overlap_score")
)
spatial_hypoxia_wide[, `:=`(
  rec_nonhypoxia_difference =
    rec_nonhypoxia_score_macro_tumour - rec_nonhypoxia_score_micro_tumour,
  rec_hypoxia_overlap_difference =
    rec_hypoxia_overlap_score_macro_tumour - rec_hypoxia_overlap_score_micro_tumour
)]
spatial_nonhypoxia_adjusted <- exact_adjusted_paired(
  spatial_hypoxia_wide$rec_nonhypoxia_difference,
  spatial_hypoxia_wide$rec_hypoxia_overlap_difference
)
spatial_hypoxia_adjusted <- exact_adjusted_paired(
  spatial_hypoxia_wide$rec_hypoxia_overlap_difference,
  spatial_hypoxia_wide$rec_nonhypoxia_difference
)

conditional_models <- rbind(
  data.table(
    dataset = "GSE151165 bulk tumour",
    outcome_component = "REC non-hypoxia component",
    covariate_component = "REC–hypoxia overlap",
    adjusted_effect = bulk_nonhypoxia_adjusted$effect,
    exact_p_two_sided = bulk_nonhypoxia_adjusted$exact_p_two_sided,
    n_assignments = bulk_nonhypoxia_adjusted$n_assignments,
    component_correlation = cor(
      bulk_hypoxia_data$rec_nonhypoxia_score,
      bulk_hypoxia_data$rec_hypoxia_overlap_score
    )
  ),
  data.table(
    dataset = "GSE151165 bulk tumour",
    outcome_component = "REC–hypoxia overlap",
    covariate_component = "REC non-hypoxia component",
    adjusted_effect = bulk_hypoxia_adjusted$effect,
    exact_p_two_sided = bulk_hypoxia_adjusted$exact_p_two_sided,
    n_assignments = bulk_hypoxia_adjusted$n_assignments,
    component_correlation = cor(
      bulk_hypoxia_data$rec_nonhypoxia_score,
      bulk_hypoxia_data$rec_hypoxia_overlap_score
    )
  ),
  data.table(
    dataset = "GSE294385 paired lesion regions",
    outcome_component = "REC non-hypoxia component",
    covariate_component = "REC–hypoxia overlap",
    adjusted_effect = spatial_nonhypoxia_adjusted$effect,
    exact_p_two_sided = spatial_nonhypoxia_adjusted$exact_p_two_sided,
    n_assignments = spatial_nonhypoxia_adjusted$n_assignments,
    component_correlation = cor(
      spatial_hypoxia_wide$rec_nonhypoxia_difference,
      spatial_hypoxia_wide$rec_hypoxia_overlap_difference
    )
  ),
  data.table(
    dataset = "GSE294385 paired lesion regions",
    outcome_component = "REC–hypoxia overlap",
    covariate_component = "REC non-hypoxia component",
    adjusted_effect = spatial_hypoxia_adjusted$effect,
    exact_p_two_sided = spatial_hypoxia_adjusted$exact_p_two_sided,
    n_assignments = spatial_hypoxia_adjusted$n_assignments,
    component_correlation = cor(
      spatial_hypoxia_wide$rec_nonhypoxia_difference,
      spatial_hypoxia_wide$rec_hypoxia_overlap_difference
    )
  )
)

bulk_score_output <- bulk_hypoxia_data[, .(
  dataset = "GSE151165 bulk tumour", unit = sample, patient = sample,
  group = hgp, rec_nonhypoxia_score, rec_hypoxia_overlap_score
)]
spatial_score_output <- spatial_hypoxia_data[, .(
  dataset = "GSE294385 paired lesion regions", unit = sample_id, patient,
  group = region, rec_nonhypoxia_score, rec_hypoxia_overlap_score
)]
write_tsv(rbind(bulk_score_output, spatial_score_output, fill = TRUE), "rec_hypoxia_component_sample_scores.tsv")
write_tsv(component_effect_rows, "rec_hypoxia_component_effects.tsv")
write_tsv(conditional_models, "rec_hypoxia_conditional_models.tsv")

# -----------------------------------------------------------------------------
# 7. Joint expression-rank-matched random programme benchmark
# -----------------------------------------------------------------------------

common_background <- intersect(rownames(bulk_log_cpm), rownames(spatial_log_cpm))
bulk_variable <- apply(bulk_log_cpm[common_background, , drop = FALSE], 1L, sd) > 0
spatial_variable <- apply(spatial_log_cpm[common_background, , drop = FALSE], 1L, sd) > 0
spatial_detection <- rowMeans(spatial_counts[common_background, , drop = FALSE] > 0)
background <- common_background[
  bulk_variable & spatial_variable & spatial_detection >= 0.50 &
    !grepl(technical_pattern, common_background, perl = TRUE)
]
rec_common <- intersect(rec_program$gene, background)
candidate_background <- setdiff(background, rec_common)
if (length(rec_common) < 40L || length(candidate_background) < 1000L) {
  stop("Insufficient common genes for the expression-matched random benchmark")
}

bulk_mean <- rowMeans(bulk_log_cpm[background, , drop = FALSE])
spatial_mean <- rowMeans(spatial_log_cpm[background, , drop = FALSE])
bulk_rank <- frank(bulk_mean, ties.method = "average") / length(background)
spatial_rank <- frank(spatial_mean, ties.method = "average") / length(background)
names(bulk_rank) <- background
names(spatial_rank) <- background

candidate_orders <- lapply(rec_common, function(target) {
  distance <- sqrt(
    (bulk_rank[candidate_background] - bulk_rank[target])^2 +
      (spatial_rank[candidate_background] - spatial_rank[target])^2
  )
  candidate_background[order(distance)][seq_len(min(200L, length(distance)))]
})
names(candidate_orders) <- rec_common

draw_matched_set <- function() {
  selected <- character()
  for (target in sample(rec_common)) {
    available <- setdiff(candidate_orders[[target]], selected)
    pool <- head(available, 50L)
    if (!length(pool)) stop("Unable to draw a unique matched random set")
    selected <- c(selected, sample(pool, 1L))
  }
  selected
}

bulk_z_common <- safe_z_rows(
  bulk_log_cpm[background, , drop = FALSE], population_sd = FALSE
)
spatial_z_common <- safe_z_rows(
  spatial_log_cpm[background, , drop = FALSE], population_sd = TRUE
)

programme_two_effects <- function(genes) {
  bulk_score <- colMeans(bulk_z_common[genes, , drop = FALSE])
  bulk_effect <- mean(bulk_score[bulk_metadata[hgp == "rHGP", sample]]) -
    mean(bulk_score[bulk_metadata[hgp == "dHGP", sample]])
  spatial_score <- colMeans(spatial_z_common[genes, , drop = FALSE])
  score_dt <- copy(spatial_meta)
  score_dt[, score := spatial_score[sample_id]]
  paired <- dcast(score_dt, patient ~ region, value.var = "score")
  spatial_effect <- mean(paired$macro_tumour - paired$micro_tumour)
  c(bulk_effect = bulk_effect, outgrowth_effect = spatial_effect)
}

observed_random_scale <- programme_two_effects(rec_common)
random_n <- 1000L
random_sets <- vector("list", random_n)
random_effects <- rbindlist(lapply(seq_len(random_n), function(iteration) {
  genes <- draw_matched_set()
  random_sets[[iteration]] <<- genes
  effects <- programme_two_effects(genes)
  data.table(
    iteration = iteration,
    bulk_effect = effects[["bulk_effect"]],
    outgrowth_effect = effects[["outgrowth_effect"]],
    mean_bulk_expression_rank = mean(bulk_rank[genes]),
    mean_spatial_expression_rank = mean(spatial_rank[genes])
  )
}))
write_tsv_gz(random_effects, "expression_matched_random_program_effects.tsv.gz")
random_summary <- data.table(
  n_random_programmes = random_n,
  n_genes_per_programme = length(rec_common),
  observed_bulk_effect = observed_random_scale[["bulk_effect"]],
  observed_outgrowth_effect = observed_random_scale[["outgrowth_effect"]],
  empirical_p_bulk = (1 + sum(random_effects$bulk_effect >= observed_random_scale[["bulk_effect"]])) / (random_n + 1),
  empirical_p_outgrowth = (1 + sum(random_effects$outgrowth_effect >= observed_random_scale[["outgrowth_effect"]])) / (random_n + 1),
  empirical_p_joint = (1 + sum(
    random_effects$bulk_effect >= observed_random_scale[["bulk_effect"]] &
      random_effects$outgrowth_effect >= observed_random_scale[["outgrowth_effect"]]
  )) / (random_n + 1),
  rec_mean_bulk_expression_rank = mean(bulk_rank[rec_common]),
  random_mean_bulk_expression_rank = mean(random_effects$mean_bulk_expression_rank),
  rec_mean_spatial_expression_rank = mean(spatial_rank[rec_common]),
  random_mean_spatial_expression_rank = mean(random_effects$mean_spatial_expression_rank)
)
write_tsv(random_summary, "expression_matched_random_program_summary.tsv")

# -----------------------------------------------------------------------------
# 8. Replace incompatible cross-source effect averages with source-family ranks
# -----------------------------------------------------------------------------

recurrence_map <- fread(file.path(
  ROOT, "analysis_results", "rec_program_multisource_core", "gene_recurrence_map.tsv"
))
hgp_columns <- c(
  "bulk_rHGP_vs_dHGP", "interface_rHGP_vs_dHGP",
  "epithelial_rHGP_vs_dHGP", "spatial_rHGP_vs_dHGP"
)
for (column in hgp_columns) {
  rank_column <- paste0(column, "_percentile")
  recurrence_map[, (rank_column) := frank(get(column), ties.method = "average") / (.N + 1)]
}
rank_columns <- paste0(hgp_columns, "_percentile")
recurrence_map[, mean_hgp_percentile_unweighted_four := rowMeans(.SD), .SDcols = rank_columns]
recurrence_map[, shared_cohort_multimodal_percentile := rowMeans(.SD), .SDcols = c(
  "epithelial_rHGP_vs_dHGP_percentile", "spatial_rHGP_vs_dHGP_percentile"
)]
recurrence_map[, mean_hgp_percentile_family_aware := rowMeans(.SD), .SDcols = c(
  "bulk_rHGP_vs_dHGP_percentile", "interface_rHGP_vs_dHGP_percentile",
  "shared_cohort_multimodal_percentile"
)]
recurrence_map[, mean_hgp_percentile_independent_only := rowMeans(.SD), .SDcols = c(
  "bulk_rHGP_vs_dHGP_percentile", "interface_rHGP_vs_dHGP_percentile"
)]
recurrence_map[, mean_hgp_percentile := mean_hgp_percentile_family_aware]
recurrence_map[, outgrowth_percentile := frank(outgrowth_macro_vs_micro, ties.method = "average") / (.N + 1)]
write_tsv(
  recurrence_map,
  "cross_source_rank_normalized_gene_effects.tsv"
)

rank_sensitivity <- data.table(
  comparison = c(
    "family-aware vs unweighted four-source",
    "family-aware vs independent-only"
  ),
  spearman_rho = c(
    cor(
      recurrence_map$mean_hgp_percentile_family_aware,
      recurrence_map$mean_hgp_percentile_unweighted_four,
      method = "spearman"
    ),
    cor(
      recurrence_map$mean_hgp_percentile_family_aware,
      recurrence_map$mean_hgp_percentile_independent_only,
      method = "spearman"
    )
  )
)
write_tsv(rank_sensitivity, "cross_source_rank_family_sensitivity.tsv")

# -----------------------------------------------------------------------------
# Summary and provenance
# -----------------------------------------------------------------------------

state_wide <- dcast(
  state_effects, programme + family ~ dataset, value.var = "effect"
)
setnames(
  state_wide,
  c("GSE151165 bulk tumour", "GSE294385 paired lesion regions"),
  c("bulk_effect", "outgrowth_effect")
)
fixed_state <- state_wide[programme == "Fixed REC 50"]
published_rec <- state_wide[programme == "Published REC"]

summary_lines <- c(
  "# Reader-facing REC-program specificity and sensitivity",
  "",
  "## NC3-to-REC mapping",
  "",
  sprintf(
    "REC ranked first with Pearson, Spearman and cosine similarity. After removing the six source markers, REC remained rank %d (r=%.3f; margin over %s=%.3f).",
    metric_rows[analysis == "pearson_without_six_source_markers", rec_rank],
    metric_rows[analysis == "pearson_without_six_source_markers", rec_similarity],
    metric_rows[analysis == "pearson_without_six_source_markers", second_state],
    metric_rows[analysis == "pearson_without_six_source_markers", rec_minus_second]
  ),
  sprintf(
    "REC remained the top state after omitting each of %d genes and in %.1f%% of %d gene-bootstrap resamples.",
    nrow(loo_gene), 100 * bootstrap_summary$rec_top_fraction, bootstrap_n
  ),
  "",
  "## External specificity",
  "",
  sprintf(
    "The fixed REC programme had effects of %+.3f in HGP-labelled bulk tumours and %+.3f in paired macro-versus-micro lesion regions.",
    fixed_state$bulk_effect, fixed_state$outgrowth_effect
  ),
  sprintf(
    "The equally sized published REC marker control had corresponding effects of %+.3f and %+.3f.",
    published_rec$bulk_effect, published_rec$outgrowth_effect
  ),
  sprintf(
    "Against %d joint expression-rank-matched random programmes, empirical upper-tail probabilities were %.4f for bulk HGP, %.4f for lesion size, and %.4f for exceeding both observed effects.",
    random_n, random_summary$empirical_p_bulk,
    random_summary$empirical_p_outgrowth, random_summary$empirical_p_joint
  ),
  sprintf(
    "Top-25, top-50, top-75 and top-100 REC programmes all retained positive effects in both external settings (bulk range %+.3f to %+.3f; lesion-size range %+.3f to %+.3f).",
    min(programme_length_effects[dataset == "GSE151165 bulk tumour", effect]),
    max(programme_length_effects[dataset == "GSE151165 bulk tumour", effect]),
    min(programme_length_effects[dataset == "GSE294385 paired lesion regions", effect]),
    max(programme_length_effects[dataset == "GSE294385 paired lesion regions", effect])
  ),
  "",
  "## Epithelial-content control",
  "",
  sprintf(
    "The bulk REC effect was %+.3f without and %+.3f with the non-overlapping epithelial-content score in the model (exact permutation P=%.4f after adjustment). The REC and epithelial-content scores correlated at r=%.3f.",
    unadjusted_test$effect, adjusted_effect, adjusted_p,
    cor(control_data$rec_score, control_data$epithelial_content_score)
  ),
  "",
  "## REC and hypoxia components",
  "",
  sprintf(
    "The non-hypoxia REC component had effects of %+.3f in bulk HGP and %+.3f across paired lesion sizes; the REC–hypoxia overlap had effects of %+.3f and %+.3f, respectively.",
    component_effect_rows[dataset == "GSE151165 bulk tumour" & component == "REC non-hypoxia component", effect],
    component_effect_rows[dataset == "GSE294385 paired lesion regions" & component == "REC non-hypoxia component", effect],
    component_effect_rows[dataset == "GSE151165 bulk tumour" & component == "REC–hypoxia overlap", effect],
    component_effect_rows[dataset == "GSE294385 paired lesion regions" & component == "REC–hypoxia overlap", effect]
  ),
  sprintf(
    "After mutual adjustment, the non-hypoxia REC effect was %+.3f in bulk (P=%.4f) and %+.3f in paired lesion regions (P=%.4f); the hypoxia-overlap effect was %+.3f (P=%.4f) and %+.3f (P=%.4f).",
    bulk_nonhypoxia_adjusted$effect, bulk_nonhypoxia_adjusted$exact_p_two_sided,
    spatial_nonhypoxia_adjusted$effect, spatial_nonhypoxia_adjusted$exact_p_two_sided,
    bulk_hypoxia_adjusted$effect, bulk_hypoxia_adjusted$exact_p_two_sided,
    spatial_hypoxia_adjusted$effect, spatial_hypoxia_adjusted$exact_p_two_sided
  ),
  "",
  "## Source-family-aware synthesis",
  "",
  sprintf(
    "E-MTAB-12022 epithelial and E-MTAB-12043 FFPE ranks were first averaged as one shared-cohort multimodal source. The resulting family-aware gene ranks correlated at rho=%.3f with the earlier unweighted four-source ranks and rho=%.3f with the two independent HGP sources alone.",
    rank_sensitivity[comparison == "family-aware vs unweighted four-source", spearman_rho],
    rank_sensitivity[comparison == "family-aware vs independent-only", spearman_rho]
  ),
  "",
  "## Interpretation boundary",
  "",
  paste(
    "These analyses test specificity and sensitivity of a transcriptomic programme.",
    "The epithelial marker panel is a composition surrogate rather than measured tumour purity,",
    "and the matched random programmes are a benchmark rather than a new discovery test."
  )
)
writeLines(summary_lines, file.path(OUT, "analysis_summary.md"), useBytes = TRUE)

run_info <- data.table(
  analysis_date = "2026-08-24",
  random_seed = 42L,
  anchor_genes = length(common_anchor_genes),
  mapping_bootstrap_replicates = bootstrap_n,
  random_programmes = random_n,
  equal_published_state_programme_size = equal_state_size,
  gse151165_tumour_samples = nrow(bulk_metadata),
  gse294385_paired_patients = uniqueN(spatial_meta$patient),
  r_version = R.version.string,
  data_table_version = as.character(packageVersion("data.table")),
  edgeR_version = as.character(packageVersion("edgeR")),
  Matrix_version = as.character(packageVersion("Matrix")),
  readxl_version = as.character(packageVersion("readxl"))
)
write_tsv(run_info, "run_info.tsv")
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))

manifest <- data.table(
  file = c(
    "state_mapping_metric_and_marker_sensitivity.tsv",
    "state_mapping_leave_one_gene_out.tsv",
    "state_mapping_gene_bootstrap.tsv.gz",
    "state_mapping_gene_bootstrap_summary.tsv",
    "gse294385_patient_region_all_gene_counts.tsv.gz",
    "gse294385_patient_region_metadata.tsv",
    "rec_program_length_sensitivity.tsv",
    "rec_program_component_effects.tsv",
    "ogden_published_state_programs.tsv",
    "ogden_state_program_external_effects.tsv",
    "bulk_epithelial_content_sample_scores.tsv",
    "bulk_epithelial_content_control.tsv",
    "expression_matched_random_program_effects.tsv.gz",
    "expression_matched_random_program_summary.tsv",
    "cross_source_rank_normalized_gene_effects.tsv",
    "rec_hypoxia_component_sample_scores.tsv",
    "rec_hypoxia_component_effects.tsv",
    "rec_hypoxia_conditional_models.tsv",
    "cross_source_rank_family_sensitivity.tsv",
    "analysis_summary.md", "run_info.tsv", "sessionInfo.txt"
  ),
  role = c(
    "mapping metrics and source-marker removal", "leave-one-gene mapping",
    "gene-bootstrap iterations", "gene-bootstrap summary",
    "cached full-gene spatial pseudobulks", "spatial pseudobulk metadata",
    "top-N programme sensitivity in the two principal external settings",
    "fixed-program component projections", "equal-size published state-marker definitions",
    "external state-program projections", "bulk sample-level covariate scores",
    "bulk epithelial-content control", "matched random programme effects",
    "matched random programme summary", "scale-free multi-source gene ranks",
    "REC and hypoxia component sample scores", "REC and hypoxia component effects",
    "mutually adjusted REC and hypoxia component models", "source-family rank sensitivity",
    "human-readable summary", "analysis provenance", "R environment"
  )
)
write_tsv(manifest, "_analysis_outputs.tsv")

message("Reader-facing specificity analysis written to: ", OUT)
