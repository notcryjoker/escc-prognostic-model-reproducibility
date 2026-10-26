# ============================================================================
# 02b_GSE53624_preprocess.R
# Purpose: preprocessing of the GSE53624 external-validation cohort (119 paired ESCC samples, Agilent GPL18109)
#   1. Parse the series matrix (238 samples = 119 tumour + 119 adjacent normal)
#   2. GPL18109 platform annotation: probes -> gene symbols (median for multiple probes)
#   3. Match clinical data (survival time + death status + TNM)
# Input: data/geo/GSE53624_series_matrix.txt.gz, GSE53624_clinical.xlsx.gz,
#       data/geo/platform_anno/GPL18109_family.soft.gz
# Output: results/geo_GSE53624.rds   tumour-sample expression (gene symbols x 119) + clinical
# ============================================================================
rm(list = ls())
suppressPackageStartupMessages({
  library(data.table); library(openxlsx); library(matrixStats)
})

base <- "/workspace/escc_bioinformatics"
out_dir <- file.path(base, "results")

# ---------- 1. Read platform annotation (probes -> gene symbols) ----------
# Note: the GEO SOFT platform table for GPL18109 contains only ID/COL/ROW/SPOT_ID/ControlType/SEQUENCE,
#       with no gene-symbol column; the Ensembl BioMart Agilent 4x44K v2 probe annotation is used instead.
bm_file   <- file.path(base, "data/geo/platform_anno/GPL18109_probe2gene_biomart.tsv")
soft_file <- file.path(base, "data/geo/platform_anno/GPL18109_family.soft.gz")

if (file.exists(bm_file)) {
  cat("Using BioMart probe annotation:", bm_file, "\n")
  probe_map <- read.delim(bm_file, stringsAsFactors = FALSE, colClasses = "character")
  names(probe_map) <- c("probe", "symbol")
  probe_map <- probe_map[!duplicated(probe_map[, c("probe", "symbol")]), ]
} else {
  cat("Fallback: parsing the SOFT platform annotation:", soft_file, "\n")
  con <- gzfile(soft_file, "rt")
  in_platform <- FALSE
  probe_map <- NULL
  repeat {
    line <- readLines(con, n = 1, warn = FALSE)
    if (length(line) == 0) break
    if (grepl("^!platform_table_begin", line)) { in_platform <- TRUE; next }
    if (grepl("^!platform_table_end", line)) break
    if (in_platform) {
      parts <- strsplit(line, "\t", fixed = TRUE)[[1]]
      if (is.null(probe_map)) {
        id_col  <- grep("^ID$", parts)
        sym_col <- grep("GeneSymbol|GENE_SYMBOL|Symbol", parts, ignore.case = TRUE)[1]
        if (is.na(sym_col)) stop("The SOFT platform table has no gene-symbol column; please provide a BioMart annotation file")
        probe_map <- data.frame(probe = character(), symbol = character())
        next
      }
      if (length(parts) >= max(id_col, sym_col)) {
        probe_map <- rbind(probe_map,
                           data.frame(probe = parts[id_col], symbol = parts[sym_col]))
      }
    }
  }
  close(con)
}
cat("Platform annotation entries:", nrow(probe_map), "\n")
probe_map <- probe_map[!is.na(probe_map$symbol) & probe_map$symbol != "" &
                         !grepl("^NA$", probe_map$symbol), ]
cat("Valid symbol entries:", nrow(probe_map), "| unique probes:", length(unique(probe_map$probe)), "\n")

# ---------- 2. Parse the series matrix ----------
sm_file <- file.path(base, "data/geo/GSE53624_series_matrix.txt.gz")
cat("Parsing series matrix:", sm_file, "\n")
lines <- readLines(gzfile(sm_file, "rt"))
# Locate the table
tbl_start <- grep("^!series_matrix_table_begin", lines)
tbl_end <- grep("^!series_matrix_table_end", lines)
header <- strsplit(lines[tbl_start + 1], "\t", fixed = TRUE)[[1]]
data_lines <- lines[(tbl_start + 2):(tbl_end - 1)]
cat("Samples:", length(header) - 1, "| data rows:", length(data_lines), "\n")

# Read in chunks into a matrix (first column is the probe ID; drop it before numeric conversion)
mat <- do.call(rbind, lapply(data_lines, function(l)
  as.numeric(strsplit(l, "\t", fixed = TRUE)[[1]])[-1]))
colnames(mat) <- header[-1]
rownames(mat) <- sapply(strsplit(data_lines, "\t", fixed = TRUE), `[`, 1)
cat("Expression matrix:", dim(mat), "\n")

# ---------- 3. Sample annotation: extract patient IDs and tissue type ----------
title_line <- grep("^!Sample_title", lines)
titles <- strsplit(sub("^!Sample_title\t", "", lines[title_line]), "\t", fixed = TRUE)[[1]]
cat("Title examples:", head(titles, 4), "\n")
# Tumour samples: title contains "cancer"; normal: contains "normal"
is_tumor <- grepl("cancer", titles, ignore.case = TRUE)
cat("Tumour samples:", sum(is_tumor), "| normal samples:", sum(!is_tumor), "\n")

# Extract patient IDs from Sample_characteristics_ch1 ("patient id: ec224")
char_lines <- grep("^!Sample_characteristics_ch1", lines, value = TRUE)
chars <- strsplit(sub("^!Sample_characteristics_ch1\t", "", char_lines[1]), "\t",
                  fixed = TRUE)[[1]]
pat <- sub(".*patient id:\\s*([A-Za-z0-9]+).*", "\\1", chars)
if (!all(grepl("^ec[0-9]+$", pat))) {
  # Fallback: extract the numeric ID from the title and prefix it with ec
  pat <- paste0("ec", sub(".*patient ([0-9]+).*", "\\1", titles))
}
cat("Patient ID examples:", head(pat, 4), "\n")

# ---------- 4. Probes -> gene symbols ----------
common <- intersect(rownames(mat), probe_map$probe)
cat("Matched probes:", length(common), "/", nrow(mat), "\n")
mat_c <- mat[common, ]
sym <- probe_map$symbol[match(common, probe_map$probe)]
expr_sym <- do.call(rbind, lapply(split(seq_len(nrow(mat_c)), sym), function(idx) {
  if (length(idx) == 1) mat_c[idx, ] else colMedians(mat_c[idx, , drop = FALSE])
}))
rownames(expr_sym) <- names(split(seq_len(nrow(mat_c)), sym))
cat("Symbol matrix:", dim(expr_sym), "\n")

# ---------- 5. Match clinical data ----------
clin <- read.xlsx(file.path(base, "data/geo/GSE53624_clinical.xlsx"))
cat("Clinical:", dim(clin), "| columns:", paste(names(clin), collapse = ", "), "\n")
tumor_cols <- is_tumor
tumor_expr <- expr_sym[, tumor_cols, drop = FALSE]
tumor_pat <- pat[tumor_cols]
colnames(tumor_expr) <- tumor_pat

# Clinical matching (Patient.ID vs IDs)
m <- match(tumor_pat, clin$Patient.ID)
clin_t <- clin[m, ]
cat("Clinical matches:", sum(!is.na(m)), "/", length(tumor_pat), "\n")

# Build the survival table: time (months -> days), status
clin_t$os_time <- as.numeric(clin_t$`Survival.time(months)`) * 30.44
clin_t$os_status <- ifelse(clin_t$`Death.at.FU` == "yes", 1, 0)
clin_t$stage <- clin_t$TNM.stage
cat("Complete survival:", sum(complete.cases(clin_t$os_time, clin_t$os_status)),
    "| events:", sum(clin_t$os_status == 1, na.rm = TRUE), "\n")
cat("Stage distribution:", paste(names(table(clin_t$stage)), table(clin_t$stage), sep = "=", collapse = ", "), "\n")

# ---------- 6. Check the 17 model genes ----------
model <- readRDS(file.path(out_dir, "lasso_model.rds"))
model_genes <- model$genes
avail_genes <- intersect(model_genes, rownames(tumor_expr))
cat("Model genes available in GSE53624:", length(avail_genes), "/", length(model_genes), "\n")
if (length(avail_genes) < length(model_genes))
  cat("Missing:", paste(setdiff(model_genes, rownames(tumor_expr)), collapse = ", "), "\n")

# ---------- 7. Save ----------
saveRDS(list(expr = tumor_expr, pheno = clin_t, probe_map_n = nrow(probe_map)),
        file.path(out_dir, "geo_GSE53624.rds"))
cat("\n02b_GSE53624_preprocess done\n")
