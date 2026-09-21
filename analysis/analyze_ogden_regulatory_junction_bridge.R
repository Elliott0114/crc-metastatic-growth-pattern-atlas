#!/usr/bin/env Rscript

# Analysis: ATAC-informed Ogden regulon-to-junction bridge
# Date: 2026-09-01
# Random seed: 42
# Key packages and versions are written to the provenance/session files.

suppressPackageStartupMessages({
  library(edgeR)
  library(readxl)
})

options(stringsAsFactors = FALSE)
set.seed(42)

root <- normalizePath(".", mustWork = TRUE)
phase2_dir <- file.path(
  root,
  "analysis_results",
  "deep_biology_upgrade_2026-08-31",
  "phase2_mechanistic_specificity"
)
out_dir <- file.path(phase2_dir, "ogden_regulatory_junction_bridge")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

spec_path <- file.path(
  root,
  "metadata",
  "ogden_atac_regulatory_junction_bridge_spec_2026-09-01.md"
)
regulon_workbook <- file.path(
  root,
  "data_sources",
  "Ogden_2025_CRLM_multiome",
  "supplementary_tables",
  "mmc8.xlsx"
)
counts_path <- file.path(phase2_dir, "ogden_direct_state_pseudobulk_counts.tsv.gz")
metadata_path <- file.path(phase2_dir, "ogden_direct_state_pseudobulk_metadata.tsv")
universe_path <- file.path(
  root,
  "analysis_results",
  "ogden_anchor_program",
  "limma_voom_all_genes.tsv.gz"
)
recurrent_path <- file.path(phase2_dir, "ogden_rec_principal_consistent_genes.tsv")
reactome_gmt <- file.path(phase2_dir, "msigdb_c2.cp.reactome_2025.1.Hs.gmt")

required <- c(
  spec_path,
  regulon_workbook,
  counts_path,
  metadata_path,
  universe_path,
  recurrent_path,
  reactome_gmt
)
if (!all(file.exists(required))) {
  stop("One or more regulatory-bridge inputs are missing")
}

write_tsv <- function(x, path) {
  write.table(
    x,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA",
    fileEncoding = "UTF-8"
  )
}

sha256 <- function(path) {
  output <- system2("sha256sum", shQuote(path), stdout = TRUE)
  strsplit(output[[1]], "[[:space:]]+")[[1]][[1]]
}

clean_genes <- function(values) {
  genes <- toupper(trimws(as.character(values)))
  unique(genes[!is.na(genes) & nzchar(genes)])
}

read_gmt_set <- function(path, set_id) {
  lines <- readLines(path, encoding = "UTF-8", warn = FALSE)
  fields <- strsplit(lines, "\t", fixed = TRUE)
  selected <- fields[vapply(fields, function(x) length(x) > 0L && x[[1]] == set_id, logical(1))]
  if (length(selected) != 1L) {
    stop(sprintf("Expected one GMT entry for %s", set_id))
  }
  clean_genes(selected[[1]][-(1:2)])
}

exact_sign_flip <- function(values) {
  values <- as.numeric(values)
  n <- length(values)
  if (n == 0L) {
    return(c(p = NA_real_, assignments = 0L))
  }
  observed <- abs(mean(values))
  assignments <- expand.grid(rep(list(c(-1, 1)), n))
  null <- apply(assignments, 1, function(signs) mean(values * signs))
  c(
    p = mean(abs(null) >= observed - 1e-15),
    assignments = nrow(assignments)
  )
}

bootstrap_mean_ci <- function(values, n_boot = 10000L) {
  draws <- replicate(n_boot, mean(sample(values, replace = TRUE)))
  unname(quantile(draws, c(0.025, 0.975), na.rm = TRUE))
}

candidate_class <- function(regulon_id) {
  if (regulon_id %in% c("AP1_COMBINED", "JUNB", "FOSL2")) {
    return("primary_AP1")
  }
  if (regulon_id %in% c("NFKB_COMBINED", "RELB")) {
    return("inflammatory_comparator")
  }
  if (regulon_id %in% c("TEAD1", "TEAD4")) {
    return("mechanotransduction_comparator")
  }
  if (regulon_id %in% c("ASCL2", "LEF1", "CDX2", "HNF4A")) {
    return("lineage_WNT_comparator")
  }
  "exploratory_source_regulon"
}

scenic <- as.data.frame(
  read_excel(regulon_workbook, sheet = "SCENIC+ TF regulons", skip = 2),
  stringsAsFactors = FALSE
)
if (ncol(scenic) != 18L) {
  stop("Unexpected number of SCENIC+ regulons")
}
scenic_sets <- lapply(scenic, clean_genes)
names(scenic_sets) <- sub("_.*$", "", names(scenic_sets))
if (anyDuplicated(names(scenic_sets))) {
  stop("SCENIC+ regulator identifiers are not unique")
}

combined <- as.data.frame(
  read_excel(regulon_workbook, sheet = "AP-1 and NFKB"),
  stringsAsFactors = FALSE
)
if (ncol(combined) != 3L) {
  stop("Unexpected combined AP-1/NF-kB table structure")
}
combined_sets <- lapply(combined, clean_genes)
names(combined_sets) <- c(
  "AP1_COMBINED",
  "AP1_EXPERIMENTAL_INTERSECTION",
  "NFKB_COMBINED"
)
regulon_sets <- c(scenic_sets, combined_sets)

regulon_long <- do.call(
  rbind,
  lapply(
    names(regulon_sets),
    function(regulon_id) {
      data.frame(
        regulon_id = regulon_id,
        candidate_class = candidate_class(regulon_id),
        gene_rank = seq_along(regulon_sets[[regulon_id]]),
        gene = regulon_sets[[regulon_id]],
        stringsAsFactors = FALSE
      )
    }
  )
)

universe_frame <- read.delim(gzfile(universe_path), check.names = FALSE)
if (!all(c("gene", "excluded_from_projection_program") %in% names(universe_frame))) {
  stop("Unexpected Ogden differential-expression universe table")
}
universe <- clean_genes(
  universe_frame$gene[!as.logical(universe_frame$excluded_from_projection_program)]
)
tight_junction <- read_gmt_set(
  reactome_gmt,
  "REACTOME_TIGHT_JUNCTION_INTERACTIONS"
)
tight_junction <- intersect(tight_junction, universe)
if (length(tight_junction) < 10L) {
  stop("Too few tight-junction genes in the tested universe")
}

ora_rows <- lapply(
  names(regulon_sets),
  function(regulon_id) {
    set_genes <- intersect(regulon_sets[[regulon_id]], universe)
    overlap <- intersect(set_genes, tight_junction)
    contingency <- matrix(
      c(
        length(overlap),
        length(tight_junction) - length(overlap),
        length(set_genes) - length(overlap),
        length(universe) - length(union(tight_junction, set_genes))
      ),
      nrow = 2,
      byrow = TRUE
    )
    fisher <- fisher.test(contingency, alternative = "greater")
    data.frame(
      regulon_id = regulon_id,
      candidate_class = candidate_class(regulon_id),
      universe_size = length(universe),
      junction_size_in_universe = length(tight_junction),
      regulon_size_in_universe = length(set_genes),
      overlap_size = length(overlap),
      odds_ratio = unname(fisher$estimate),
      fisher_one_sided_p = fisher$p.value,
      overlap_genes = paste(sort(overlap), collapse = ";"),
      stringsAsFactors = FALSE
    )
  }
)
ora <- do.call(rbind, ora_rows)
ora$fdr_bh_all_regulons <- p.adjust(ora$fisher_one_sided_p, method = "BH")
ora <- ora[order(ora$fdr_bh_all_regulons, -ora$overlap_size), , drop = FALSE]

recurrent <- read.delim(recurrent_path, check.names = FALSE)
recurrent_up <- clean_genes(recurrent$gene[recurrent$direction == "REC_up_all_three"])
triple_rows <- lapply(
  names(regulon_sets),
  function(regulon_id) {
    genes <- sort(intersect(intersect(regulon_sets[[regulon_id]], tight_junction), recurrent_up))
    if (length(genes) == 0L) {
      return(NULL)
    }
    data.frame(
      regulon_id = regulon_id,
      candidate_class = candidate_class(regulon_id),
      gene = genes,
      stringsAsFactors = FALSE
    )
  }
)
triple_bridge <- do.call(rbind, triple_rows[!vapply(triple_rows, is.null, logical(1))])
if (is.null(triple_bridge)) {
  triple_bridge <- data.frame(
    regulon_id = character(),
    candidate_class = character(),
    gene = character()
  )
}

counts_frame <- read.delim(gzfile(counts_path), check.names = FALSE)
if (!"gene" %in% names(counts_frame)) {
  stop("Pseudobulk counts table lacks gene")
}
genes <- clean_genes(counts_frame$gene)
if (length(genes) != nrow(counts_frame)) {
  stop("Pseudobulk gene identifiers are missing or duplicated")
}
counts <- as.matrix(counts_frame[, setdiff(names(counts_frame), "gene"), drop = FALSE])
storage.mode(counts) <- "double"
rownames(counts) <- toupper(counts_frame$gene)
metadata <- read.delim(metadata_path, check.names = FALSE)
if (!setequal(colnames(counts), metadata$sample_id)) {
  stop("Pseudobulk counts and metadata identifiers differ")
}
metadata <- metadata[match(colnames(counts), metadata$sample_id), , drop = FALSE]
eligibility_text <- tolower(trimws(as.character(metadata$eligible_min_cells)))
if (!all(eligibility_text %in% c("true", "false"))) {
  stop("eligible_min_cells contains values other than True/False")
}
metadata$eligible_min_cells <- eligibility_text == "true"
selected <- metadata$eligible_min_cells & metadata$state %in% c("REC", "Hypoxia", "UPR", "iREC")
selected_metadata <- metadata[selected, , drop = FALSE]
selected_counts <- counts[, selected, drop = FALSE]
log_cpm <- cpm(selected_counts, log = TRUE, prior.count = 0.5)
gene_sd <- apply(log_cpm, 1, sd)
valid_gene <- is.finite(gene_sd) & gene_sd > 0
gene_z <- sweep(log_cpm[valid_gene, , drop = FALSE], 1, rowMeans(log_cpm[valid_gene, , drop = FALSE]), "-")
gene_z <- sweep(gene_z, 1, gene_sd[valid_gene], "/")

score_rows <- list()
score_index <- 1L
for (regulon_id in names(regulon_sets)) {
  available <- intersect(regulon_sets[[regulon_id]], rownames(gene_z))
  if (length(available) < 5L) {
    next
  }
  score <- colMeans(gene_z[available, , drop = FALSE])
  score_rows[[score_index]] <- data.frame(
    selected_metadata,
    regulon_id = regulon_id,
    candidate_class = candidate_class(regulon_id),
    detected_target_genes = length(available),
    regulon_gene_z_score = as.numeric(score),
    stringsAsFactors = FALSE
  )
  score_index <- score_index + 1L
}
scores <- do.call(rbind, score_rows)

effect_rows <- list()
patient_rows <- list()
effect_index <- 1L
patient_index <- 1L
for (regulon_id in unique(scores$regulon_id)) {
  regulator <- scores[scores$regulon_id == regulon_id, , drop = FALSE]
  rec <- regulator[regulator$state == "REC", c("patient", "regulon_gene_z_score")]
  names(rec)[[2]] <- "rec_score"
  for (comparator in c("Hypoxia", "UPR", "iREC")) {
    other <- regulator[
      regulator$state == comparator,
      c("patient", "regulon_gene_z_score")
    ]
    names(other)[[2]] <- "comparator_score"
    paired <- merge(rec, other, by = "patient")
    if (nrow(paired) < 5L) {
      next
    }
    paired$difference <- paired$rec_score - paired$comparator_score
    for (row in seq_len(nrow(paired))) {
      patient_rows[[patient_index]] <- data.frame(
        regulon_id = regulon_id,
        candidate_class = candidate_class(regulon_id),
        comparator_state = comparator,
        patient = paired$patient[[row]],
        rec_score = paired$rec_score[[row]],
        comparator_score = paired$comparator_score[[row]],
        rec_minus_comparator = paired$difference[[row]],
        stringsAsFactors = FALSE
      )
      patient_index <- patient_index + 1L
    }
    ci <- bootstrap_mean_ci(paired$difference)
    exact <- exact_sign_flip(paired$difference)
    effect_rows[[effect_index]] <- data.frame(
      regulon_id = regulon_id,
      candidate_class = candidate_class(regulon_id),
      comparator_state = comparator,
      n_paired_patients = nrow(paired),
      n_positive = sum(paired$difference > 0),
      n_negative = sum(paired$difference < 0),
      mean_rec_minus_comparator = mean(paired$difference),
      median_rec_minus_comparator = median(paired$difference),
      bootstrap_ci_lower = ci[[1]],
      bootstrap_ci_upper = ci[[2]],
      exact_sign_flip_p = exact[["p"]],
      n_exact_assignments = exact[["assignments"]],
      stringsAsFactors = FALSE
    )
    effect_index <- effect_index + 1L
  }
}
effects <- do.call(rbind, effect_rows)
effects$fdr_bh_within_comparator <- ave(
  effects$exact_sign_flip_p,
  effects$comparator_state,
  FUN = function(x) p.adjust(x, method = "BH")
)
effects <- effects[order(effects$comparator_state, -effects$mean_rec_minus_comparator), , drop = FALSE]
patient_effects <- do.call(rbind, patient_rows)

write_tsv(regulon_long, file.path(out_dir, "ogden_source_regulons_long.tsv"))
write_tsv(ora, file.path(out_dir, "ogden_regulon_tight_junction_ora.tsv"))
write_tsv(triple_bridge, file.path(out_dir, "ogden_recurrent_rec_junction_regulon_bridge.tsv"))
write_tsv(scores, file.path(out_dir, "ogden_regulon_patient_state_scores.tsv"))
write_tsv(patient_effects, file.path(out_dir, "ogden_regulon_paired_patient_effects.tsv"))
write_tsv(effects, file.path(out_dir, "ogden_regulon_paired_state_summary.tsv"))

provenance <- data.frame(
  field = c(
    "analysis_date",
    "seed",
    "organism",
    "gene_namespace",
    "ora_background",
    "pseudobulk_score",
    "bootstrap_replicates",
    "specification_sha256",
    "regulon_workbook_sha256",
    "counts_sha256",
    "metadata_sha256",
    "universe_sha256",
    "recurrent_results_sha256",
    "reactome_gmt_sha256"
  ),
  value = c(
    "2026-09-01",
    "42",
    "Homo sapiens",
    "HGNC gene symbols",
    sprintf("%d tested non-technical genes", length(universe)),
    "mean gene-wise z score of log2 CPM",
    "10000",
    sha256(spec_path),
    sha256(regulon_workbook),
    sha256(counts_path),
    sha256(metadata_path),
    sha256(universe_path),
    sha256(recurrent_path),
    sha256(reactome_gmt)
  ),
  stringsAsFactors = FALSE
)
write_tsv(provenance, file.path(out_dir, "ogden_regulatory_bridge_provenance.tsv"))
sink(file.path(out_dir, "ogden_regulatory_bridge_sessionInfo.txt"))
sessionInfo()
sink()

manifest <- c(
  "# Analysis Outputs",
  "",
  "Generated: 2026-09-01  ",
  "Study type: source-constrained regulon ORA and patient-paired state scoring",
  "",
  "## Tables",
  "",
  "- `ogden_source_regulons_long.tsv` -- Source SCENIC+ and combined regulon membership.",
  "- `ogden_regulon_tight_junction_ora.tsv` -- Junction-set ORA with tested-gene background.",
  "- `ogden_recurrent_rec_junction_regulon_bridge.tsv` -- Recurrent REC-up junction targets by regulon.",
  "- `ogden_regulon_patient_state_scores.tsv` -- Patient-state regulon scores.",
  "- `ogden_regulon_paired_patient_effects.tsv` -- Patient-paired REC contrasts.",
  "- `ogden_regulon_paired_state_summary.tsv` -- Bootstrap intervals and exact sign-flip tests.",
  "",
  "## Provenance",
  "",
  "- `ogden_regulatory_bridge_provenance.tsv` and `ogden_regulatory_bridge_sessionInfo.txt`."
)
writeLines(manifest, file.path(out_dir, "_analysis_outputs.md"), useBytes = TRUE)

cat("\nRegulon-to-tight-junction ORA\n")
print(ora[, c(
  "regulon_id",
  "candidate_class",
  "regulon_size_in_universe",
  "overlap_size",
  "odds_ratio",
  "fdr_bh_all_regulons",
  "overlap_genes"
)], row.names = FALSE)
cat("\nPatient-paired candidate effects\n")
print(
  effects[effects$candidate_class != "exploratory_source_regulon", ],
  row.names = FALSE
)
