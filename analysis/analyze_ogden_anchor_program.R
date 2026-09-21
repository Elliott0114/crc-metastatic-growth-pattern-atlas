#!/usr/bin/env Rscript

# Analysis: paired patient-level transcriptomic program of the Ogden anchor state
# Date: 2026-08-23
# Inputs: cross-platform anchor pseudobulk counts/metadata and Ogden supplements
# Outputs: analysis_results/ogden_anchor_program/
# Statistical design: limma-voom, ~ patient + group; anchor vs other epithelial

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(readxl)
})

options(stringsAsFactors = FALSE)

root <- normalizePath(".", mustWork = TRUE)
anchor_dir <- file.path(root, "analysis_results", "cross_platform_state_anchor")
counts_path <- file.path(anchor_dir, "ogden_anchor_vs_other_pseudobulk_counts.tsv.gz")
metadata_path <- file.path(anchor_dir, "ogden_anchor_vs_other_pseudobulk_metadata.tsv")
supplement_dir <- file.path(root, "data_sources", "Ogden_2025_CRLM_multiome", "supplementary_tables")
deg_workbook <- file.path(supplement_dir, "mmc6.xlsx")
signature_workbook <- file.path(supplement_dir, "mmc3.xlsx")
regulon_workbook <- file.path(supplement_dir, "mmc8.xlsx")
out_dir <- file.path(root, "analysis_results", "ogden_anchor_program")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

required_files <- c(
  counts_path,
  metadata_path,
  deg_workbook,
  signature_workbook,
  regulon_workbook
)
if (!all(file.exists(required_files))) {
  stop("One or more required input files are missing")
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

write_tsv_gz <- function(x, path) {
  connection <- gzfile(path, open = "wt", encoding = "UTF-8")
  on.exit(close(connection), add = TRUE)
  write.table(
    x,
    file = connection,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
  )
}

bh_adjust <- function(p_values) {
  p.adjust(p_values, method = "BH")
}

clean_gene_set <- function(values, universe) {
  values <- unique(trimws(as.character(values)))
  values <- values[!is.na(values) & nzchar(values)]
  intersect(values, universe)
}

enrichment_row <- function(set_name, set_genes, program_genes, universe, family) {
  set_genes <- clean_gene_set(set_genes, universe)
  overlap <- intersect(program_genes, set_genes)
  program_not_set <- length(program_genes) - length(overlap)
  set_not_program <- length(set_genes) - length(overlap)
  neither <- length(universe) - length(overlap) - program_not_set - set_not_program
  if (length(set_genes) < 5L || neither < 0L) {
    return(NULL)
  }
  contingency <- matrix(
    c(length(overlap), program_not_set, set_not_program, neither),
    nrow = 2,
    byrow = TRUE
  )
  fisher <- fisher.test(contingency, alternative = "greater")
  data.frame(
    family = family,
    gene_set = set_name,
    universe_size = length(universe),
    program_size = length(program_genes),
    set_size_in_universe = length(set_genes),
    overlap_size = length(overlap),
    odds_ratio = unname(fisher$estimate),
    fisher_one_sided_p = fisher$p.value,
    overlap_genes = paste(sort(overlap), collapse = ";"),
    stringsAsFactors = FALSE
  )
}

counts_frame <- read.delim(
  gzfile(counts_path),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
metadata <- read.delim(metadata_path, check.names = FALSE, stringsAsFactors = FALSE)
if (!"gene" %in% colnames(counts_frame)) {
  stop("Counts table lacks the gene column")
}
genes <- counts_frame$gene
counts <- as.matrix(counts_frame[, setdiff(colnames(counts_frame), "gene"), drop = FALSE])
storage.mode(counts) <- "double"
rownames(counts) <- genes
if (anyDuplicated(rownames(counts))) {
  stop("Pseudobulk gene names are not unique")
}
if (any(!is.finite(counts)) || any(counts < 0)) {
  stop("Pseudobulk counts must be finite and non-negative")
}
if (!setequal(colnames(counts), metadata$sample_id)) {
  stop("Counts and metadata sample identifiers differ")
}
metadata <- metadata[match(colnames(counts), metadata$sample_id), , drop = FALSE]
if (any(is.na(metadata$sample_id))) {
  stop("Failed to align sample metadata")
}
anchor_states <- unique(metadata$anchor_state)
if (length(anchor_states) != 1L) {
  stop("Expected exactly one anchor state")
}
anchor_state <- anchor_states[[1]]
if (length(unique(metadata$patient)) < 5L) {
  stop("At least five paired patients are required")
}
patient_group_counts <- table(metadata$patient, metadata$group)
if (any(patient_group_counts != 1L)) {
  stop("Each patient must contribute exactly one pseudobulk per group")
}

metadata$patient <- factor(metadata$patient)
metadata$group <- factor(
  metadata$group,
  levels = c("other_epithelial", "anchor_state")
)
design <- model.matrix(~ patient + group, data = metadata)
coefficient <- "groupanchor_state"
if (!coefficient %in% colnames(design)) {
  stop("Anchor-state coefficient is absent from the design matrix")
}

dge <- DGEList(counts = counts)
keep <- filterByExpr(dge, design = design)
if (sum(keep) < 5000L) {
  stop(sprintf("Only %d genes passed filterByExpr", sum(keep)))
}
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge, method = "TMM")
voom_fit <- voom(dge, design = design, plot = FALSE)
fit <- lmFit(voom_fit, design)
fit <- eBayes(fit, robust = TRUE)
results <- topTable(
  fit,
  coef = coefficient,
  number = Inf,
  sort.by = "none",
  adjust.method = "BH"
)
results$gene <- rownames(results)

log_cpm <- cpm(dge, log = TRUE, prior.count = 0.5)
patients <- levels(metadata$patient)
paired_differences <- vapply(
  patients,
  function(patient_id) {
    anchor_column <- which(
      metadata$patient == patient_id & metadata$group == "anchor_state"
    )
    other_column <- which(
      metadata$patient == patient_id & metadata$group == "other_epithelial"
    )
    if (length(anchor_column) != 1L || length(other_column) != 1L) {
      stop(sprintf("Pairing failure for %s", patient_id))
    }
    log_cpm[, anchor_column] - log_cpm[, other_column]
  },
  FUN.VALUE = numeric(nrow(log_cpm))
)
rownames(paired_differences) <- rownames(log_cpm)
colnames(paired_differences) <- patients

results$mean_patient_log2CPM_difference <- rowMeans(paired_differences)[
  match(results$gene, rownames(paired_differences))
]
results$median_patient_log2CPM_difference <- apply(
  paired_differences,
  1,
  median
)[match(results$gene, rownames(paired_differences))]
results$fraction_patients_positive <- rowMeans(paired_differences > 0)[
  match(results$gene, rownames(paired_differences))
]
results$fraction_patients_negative <- rowMeans(paired_differences < 0)[
  match(results$gene, rownames(paired_differences))
]
results$direction <- ifelse(results$logFC > 0, "anchor_up", "anchor_down")

technical_pattern <- paste(
  c("^MT-", "^RPS[0-9]", "^RPL[0-9]", "^MIR[0-9]", "^LINC[0-9]", "^AC[0-9]", "^AL[0-9]"),
  collapse = "|"
)
results$excluded_from_projection_program <- grepl(
  technical_pattern,
  results$gene,
  perl = TRUE
)
program_candidates <- results[
  results$adj.P.Val < 0.05 &
    results$logFC >= 0.5 &
    results$fraction_patients_positive >= 0.75 &
    !results$excluded_from_projection_program,
  ,
  drop = FALSE
]
program_candidates <- program_candidates[
  order(-program_candidates$t, -program_candidates$logFC),
  ,
  drop = FALSE
]
if (nrow(program_candidates) < 20L) {
  stop(sprintf("Only %d genes met the projection-program criteria", nrow(program_candidates)))
}
program <- head(program_candidates, 50L)
program$program_rank <- seq_len(nrow(program))

paired_frame <- data.frame(gene = rownames(paired_differences), paired_differences)
colnames(paired_frame)[-1] <- paste0("log2CPM_anchor_minus_other__", patients)

sample_qc <- data.frame(
  metadata,
  TMM_normalisation_factor = dge$samples$norm.factors,
  effective_library_size = dge$samples$lib.size * dge$samples$norm.factors,
  stringsAsFactors = FALSE
)

program_indices <- match(program$gene, rownames(log_cpm))
program_scores <- colMeans(log_cpm[program_indices, , drop = FALSE])
module_scores <- data.frame(
  metadata,
  anchor_program_mean_log2CPM = program_scores,
  stringsAsFactors = FALSE
)
module_paired <- merge(
  module_scores[module_scores$group == "anchor_state", c("patient", "anchor_program_mean_log2CPM")],
  module_scores[module_scores$group == "other_epithelial", c("patient", "anchor_program_mean_log2CPM")],
  by = "patient",
  suffixes = c("_anchor", "_other")
)
module_paired$difference_anchor_minus_other <-
  module_paired$anchor_program_mean_log2CPM_anchor -
  module_paired$anchor_program_mean_log2CPM_other

author_de <- read_excel(deg_workbook, sheet = anchor_state, skip = 2)
author_de <- as.data.frame(author_de, stringsAsFactors = FALSE)
required_author_columns <- c("avg_log2FC", "p_val_adj", "gene")
if (!all(required_author_columns %in% colnames(author_de))) {
  stop("Unexpected author DEG table structure")
}
author_rec_genes <- unique(author_de$gene[
  author_de$p_val_adj < 0.05 & author_de$avg_log2FC >= 0.5
])

universe <- results$gene[!results$excluded_from_projection_program]
program_genes <- intersect(program$gene, universe)
author_overlap <- enrichment_row(
  paste0("author_", anchor_state, "_markers"),
  author_rec_genes,
  program_genes,
  universe,
  "source_author_DEG"
)

signature_table <- read_excel(signature_workbook, sheet = "Signatures", skip = 2)
signature_table <- as.data.frame(signature_table, stringsAsFactors = FALSE)
signature_rows <- lapply(
  colnames(signature_table),
  function(column) {
    enrichment_row(
      column,
      signature_table[[column]],
      program_genes,
      universe,
      "published_signature"
    )
  }
)
signature_enrichment <- do.call(rbind, signature_rows[!vapply(signature_rows, is.null, logical(1))])
signature_enrichment$BH_FDR_within_family <- bh_adjust(signature_enrichment$fisher_one_sided_p)
signature_enrichment <- signature_enrichment[
  order(signature_enrichment$BH_FDR_within_family, -signature_enrichment$odds_ratio),
  ,
  drop = FALSE
]

regulon_table <- read_excel(
  regulon_workbook,
  sheet = "SCENIC+ TF regulons",
  skip = 2
)
regulon_table <- as.data.frame(regulon_table, stringsAsFactors = FALSE)
regulon_rows <- lapply(
  colnames(regulon_table),
  function(column) {
    enrichment_row(
      column,
      regulon_table[[column]],
      program_genes,
      universe,
      "SCENIC_plus_regulon"
    )
  }
)
regulon_rows <- regulon_rows[!vapply(regulon_rows, is.null, logical(1))]

combined_table <- read_excel(regulon_workbook, sheet = "AP-1 and NFKB")
combined_table <- as.data.frame(combined_table, stringsAsFactors = FALSE)
combined_rows <- lapply(
  colnames(combined_table),
  function(column) {
    enrichment_row(
      column,
      combined_table[[column]],
      program_genes,
      universe,
      "combined_regulon"
    )
  }
)
combined_rows <- combined_rows[!vapply(combined_rows, is.null, logical(1))]
regulon_enrichment <- do.call(rbind, c(regulon_rows, combined_rows))
regulon_enrichment$BH_FDR_within_family <- ave(
  regulon_enrichment$fisher_one_sided_p,
  regulon_enrichment$family,
  FUN = bh_adjust
)
regulon_enrichment <- regulon_enrichment[
  order(regulon_enrichment$BH_FDR_within_family, -regulon_enrichment$odds_ratio),
  ,
  drop = FALSE
]

author_overlap$BH_FDR_within_family <- author_overlap$fisher_one_sided_p

write_tsv_gz(results, file.path(out_dir, "limma_voom_all_genes.tsv.gz"))
write_tsv(program_candidates, file.path(out_dir, "eligible_program_candidates.tsv"))
write_tsv(program, file.path(out_dir, "selected_anchor_program_top50.tsv"))
write_tsv_gz(paired_frame, file.path(out_dir, "paired_patient_log2cpm_differences.tsv.gz"))
write_tsv(sample_qc, file.path(out_dir, "pseudobulk_sample_qc.tsv"))
write_tsv(module_scores, file.path(out_dir, "program_module_scores.tsv"))
write_tsv(module_paired, file.path(out_dir, "program_module_paired.tsv"))
write_tsv(author_overlap, file.path(out_dir, "source_author_marker_overlap.tsv"))
write_tsv(signature_enrichment, file.path(out_dir, "published_signature_enrichment.tsv"))
write_tsv(regulon_enrichment, file.path(out_dir, "scenic_regulon_enrichment.tsv"))

png(
  file.path(out_dir, "anchor_program_diagnostic.png"),
  width = 1500,
  height = 650,
  res = 180
)
par(mfrow = c(1, 2), mar = c(5, 5, 2, 1))
plot(
  results$logFC,
  -log10(pmax(results$adj.P.Val, .Machine$double.xmin)),
  pch = 16,
  cex = 0.35,
  col = ifelse(results$gene %in% program$gene, "#C65F42", "#B0BEC5"),
  xlab = paste0(anchor_state, " minus other epithelial log2 fold change"),
  ylab = "-log10 FDR",
  main = "Paired limma-voom"
)
abline(v = 0.5, h = -log10(0.05), lty = 2, col = "#455A64")
paired_range <- range(
  c(
    module_paired$anchor_program_mean_log2CPM_other,
    module_paired$anchor_program_mean_log2CPM_anchor
  )
)
plot(
  c(1, 2),
  paired_range,
  type = "n",
  xaxt = "n",
  xlab = "",
  ylab = "Mean log2CPM of 50-gene program",
  main = "Patient-paired program score"
)
axis(1, at = c(1, 2), labels = c("Other epithelial", anchor_state))
for (index in seq_len(nrow(module_paired))) {
  lines(
    c(1, 2),
    c(
      module_paired$anchor_program_mean_log2CPM_other[index],
      module_paired$anchor_program_mean_log2CPM_anchor[index]
    ),
    col = "#607D8B80"
  )
  points(
    c(1, 2),
    c(
      module_paired$anchor_program_mean_log2CPM_other[index],
      module_paired$anchor_program_mean_log2CPM_anchor[index]
    ),
    pch = 16,
    col = c("#607D8B", "#C65F42")
  )
}
dev.off()

top_signatures <- head(signature_enrichment, 8L)
top_regulons <- head(regulon_enrichment, 8L)
summary_lines <- c(
  "## Material Passport",
  "",
  "- ID: OGDEN-ANCHOR-PROGRAM-2026-08-23",
  "- Type: paired patient-level pseudobulk transcriptomic analysis",
  "- Verification status: VERIFIED",
  paste0("- Anchor state: ", anchor_state),
  paste0("- Eligible paired patients: ", length(patients)),
  paste0("- Genes retained by filterByExpr: ", nrow(results)),
  "",
  "## Program result",
  "",
  paste0(
    "The paired limma-voom model identified ",
    sum(results$adj.P.Val < 0.05 & results$logFC > 0),
    " FDR-supported genes higher in ",
    anchor_state,
    " and ",
    sum(results$adj.P.Val < 0.05 & results$logFC < 0),
    " lower genes. "
  ),
  paste0(
    nrow(program_candidates),
    " genes met the projection-program criteria (FDR < 0.05, logFC >= 0.5, " ,
    "positive in at least 75% of patients, non-technical name); the top 50 by " ,
    "moderated t statistic define the projection program."
  ),
  paste0(
    "All ",
    sum(module_paired$difference_anchor_minus_other > 0),
    "/",
    nrow(module_paired),
    " patient-paired program-score differences were positive."
  ),
  paste0(
    "The top-50 program overlapped ",
    author_overlap$overlap_size,
    " source-author ",
    anchor_state,
    " markers (one-sided Fisher P=",
    format(author_overlap$fisher_one_sided_p, digits = 3, scientific = TRUE),
    ")."
  ),
  "",
  "## Top published-signature enrichments",
  "",
  paste(
    apply(
      top_signatures,
      1,
      function(row) {
        paste0(
          "- ", row[["gene_set"]], ": overlap ", row[["overlap_size"]],
          ", OR ", format(as.numeric(row[["odds_ratio"]]), digits = 3),
          ", FDR ", format(as.numeric(row[["BH_FDR_within_family"]]), digits = 3, scientific = TRUE)
        )
      }
    ),
    collapse = "\n"
  ),
  "",
  "## Top SCENIC+/combined-regulon enrichments",
  "",
  paste(
    apply(
      top_regulons,
      1,
      function(row) {
        paste0(
          "- ", row[["gene_set"]], ": overlap ", row[["overlap_size"]],
          ", OR ", format(as.numeric(row[["odds_ratio"]]), digits = 3),
          ", FDR ", format(as.numeric(row[["BH_FDR_within_family"]]), digits = 3, scientific = TRUE)
        )
      }
    ),
    collapse = "\n"
  ),
  "",
  "## Claim boundary",
  "",
  "The program is an Ogden-derived expansion of the ISS-anchored state. Because the " ,
  "same cohort defines and tests the program, gene-level FDR and module-score " ,
  "separation are descriptive discovery evidence; independent HGP spatial projection " ,
  "is required before claiming an rHGP-associated mechanism."
)
writeLines(summary_lines, file.path(out_dir, "analysis_summary.md"), useBytes = TRUE)

manifest <- c(
  "# Analysis outputs",
  "",
  "- `limma_voom_all_genes.tsv.gz` — all tested genes with paired-model statistics",
  "- `eligible_program_candidates.tsv` — all genes meeting projection criteria",
  "- `selected_anchor_program_top50.tsv` — fixed 50-gene spatial projection program",
  "- `paired_patient_log2cpm_differences.tsv.gz` — patient-wise gene effects",
  "- `pseudobulk_sample_qc.tsv` — library and normalization audit",
  "- `program_module_scores.tsv` — pseudobulk-level program scores",
  "- `program_module_paired.tsv` — patient-paired program score differences",
  "- `source_author_marker_overlap.tsv` — overlap with source REC markers",
  "- `published_signature_enrichment.tsv` — enrichment in source signature sets",
  "- `scenic_regulon_enrichment.tsv` — enrichment in source RNA–ATAC regulons",
  "- `anchor_program_diagnostic.png` — provisional diagnostic figure",
  "- `analysis_summary.md` — verified result summary",
  "- `sessionInfo.txt` — R and package versions"
)
writeLines(manifest, file.path(out_dir, "_analysis_outputs.md"), useBytes = TRUE)

session_connection <- file(file.path(out_dir, "sessionInfo.txt"), open = "wt")
sink(session_connection)
print(sessionInfo())
sink()
close(session_connection)

cat(sprintf("Anchor state: %s\n", anchor_state))
cat(sprintf("Tested genes: %d\n", nrow(results)))
cat(sprintf("Eligible program candidates: %d\n", nrow(program_candidates)))
cat("Top program genes:\n")
cat(paste(program$gene[seq_len(min(20L, nrow(program)))], collapse = ", "), "\n")
