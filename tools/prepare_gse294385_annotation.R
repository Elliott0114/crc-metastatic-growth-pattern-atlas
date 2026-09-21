#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript --vanilla tools/prepare_gse294385_annotation.R SOURCE.rds INPUTS")
}
expected <- "7c6d6d3fe32daf0136176fec53aef180e147b4c4372bd3dd5847adb83181c609"
actual <- strsplit(system2("sha256sum", shQuote(args[1]), stdout = TRUE), " ")[[1]][1]
if (!identical(actual, expected)) stop("Source annotation checksum mismatch")
all_args <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", all_args[grep("^--file=", all_args)][1])
root <- normalizePath(file.path(dirname(script), ".."))
manifest <- read.delim(file.path(root, "data/inputs/metadata/gse294385_liver_sample_manifest.tsv"))
annotation <- readRDS(args[1])
annotation$sample <- sub("~.*$", "", annotation$Barcode)
annotation$spot_barcode <- sub("^[^~]+~", "", annotation$Barcode)
annotation <- annotation[annotation$Layer1 == "Liver", , drop = FALSE]
annotation$patient <- manifest$patient[match(annotation$sample, manifest$sample)]
if (anyNA(annotation$patient)) stop("Unmapped liver sample in source annotation")
output <- file.path(args[2], "data_sources/Liu_2026_GSE294385/visium_liver_meta_after_qc.tsv.gz")
if (file.exists(output)) stop("Output already exists; use a new input directory")
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
connection <- gzfile(output, open = "wt", encoding = "UTF-8")
write.table(annotation[, c("sample", "patient", "spot_barcode", "Barcode", "Layer1", "Layer2", "Layer3")],
            connection, sep = "\t", row.names = FALSE, quote = FALSE)
close(connection)
expected_output <- "cfcab046fe13bbc26cfb2646e1c646a28cd19493219b3746f64cefe7b52aa1bc"
actual_output <- strsplit(system2("sha256sum", shQuote(output), stdout = TRUE), " ")[[1]][1]
if (!identical(actual_output, expected_output)) stop("Prepared annotation checksum mismatch")
message("Prepared and verified ", nrow(annotation), " liver spots: ", output)
