# ============================================================================
# 02_GEO_preprocess.R
# Purpose: parse the GEO ESCC cohorts (GSE20347/GSE45670/GSE75241/GSE53624),
#       complete probe->gene-symbol mapping, multi-probe collapsing, phenotype extraction (with survival); save RDS
# Output: results/geo_<GSE>.rds  one list(expr, pheno, meta) per cohort
# ============================================================================
rm(list = ls())
suppressPackageStartupMessages({
  library(data.table); library(matrixStats); library(jsonlite)
})

base <- "/workspace/escc_bioinformatics"
geo_dir <- file.path(base, "data/geo")
anno_dir <- file.path(geo_dir, "platform_anno")
out_dir <- file.path(base, "results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---------- Utility functions ----------
read_series_matrix <- function(gzfile) {
  lines <- readLines(gzfile, warn = FALSE)
  # Locate the data table
  begin <- grep("^!series_matrix_table_begin", lines)
  end   <- grep("^!series_matrix_table_end", lines)
  header <- grep("^!Sample_", lines)
  # Extract sample-level metadata (!Sample_ rows)
  sample_meta <- list()
  for (ln in header) {
    key <- sub("^!Sample_([a-zA-Z0-9_]+).*", "\\1", lines[ln])
    vals <- strsplit(sub("^![^\t]+\t", "", lines[ln]), "\t", fixed = TRUE)[[1]]
    vals <- gsub('^"|"$', "", vals)
    sample_meta[[key]] <- vals
  }
  # Read the data table (GEO series-matrix fields are double-quoted; parse with the default quote='"')
  tab <- read.table(text = paste(lines[(begin + 1):(end - 1)], collapse = "\n"),
                    sep = "\t", header = TRUE, check.names = FALSE,
                    comment.char = "", stringsAsFactors = FALSE)
  list(expr = tab, meta = sample_meta, platform = unique(sample_meta$platform_id))
}

collapse_probes <- function(expr_mat, probe2gene) {
  # expr_mat: probes x samples; probe2gene: data.frame(probe, gene)
  expr_mat <- as.matrix(expr_mat)
  mode(expr_mat) <- "numeric"
  m <- merge(expr_mat, probe2gene, by.x = 0, by.y = "probe", all.x = FALSE)
  genes <- unique(m$gene)
  out <- do.call(rbind, lapply(genes, function(g) {
    rows <- m[m$gene == g, , drop = FALSE]
    vals <- as.matrix(rows[, colnames(expr_mat), drop = FALSE])
    mode(vals) <- "numeric"
    if (nrow(rows) == 1) as.numeric(vals[1, ]) else colMeans(vals, na.rm = TRUE)
  }))
  rownames(out) <- genes
  colnames(out) <- colnames(expr_mat)
  out
}

parse_annot <- function(path, type = c("annot", "soft")) {
  type <- match.arg(type)
  con <- gzfile(path, "rt")
  on.exit(close(con))
  if (type == "annot") {
    lines <- readLines(con)
    b <- grep("^!platform_table_begin", lines); e <- grep("^!platform_table_end", lines)
    tab <- read.table(text = paste(lines[(b + 1):(e - 1)], collapse = "\n"),
                      sep = "\t", header = TRUE, quote = "", comment.char = "",
                      stringsAsFactors = FALSE)
    probe_col <- intersect(c("ID", "SPOT_ID"), colnames(tab))[1]
    gene_col <- grep("^Gene Symbol$|^GENE_SYMBOL$|^Symbol$", colnames(tab), value = TRUE)[1]
    if (is.na(gene_col)) gene_col <- grep("ymbol", colnames(tab), value = TRUE)[1]
    out <- data.frame(probe = tab[[probe_col]], gene = trimws(tab[[gene_col]]),
                      stringsAsFactors = FALSE)
  } else {
    # family.soft.gz: read the data table near ^#ID_REF
    lines <- readLines(con)
    idx <- grep("^!platform_table_begin", lines)
    if (length(idx) == 0) idx <- grep("^#ID_REF", lines) - 1
    e <- grep("^!platform_table_end", lines)
    if (length(e) == 0) e <- idx + 1 + which(lines[(idx + 2):length(lines)] == "")
    tab <- read.table(text = paste(lines[(idx + 1):(min(e[1], length(lines)) - 1)], collapse = "\n"),
                      sep = "\t", header = TRUE, quote = "", comment.char = "",
                      stringsAsFactors = FALSE)
    probe_col <- grep("^ID$|^ID_REF$", colnames(tab), value = TRUE)[1]
    gene_col <- grep("ymbol|GENE", colnames(tab), value = TRUE)
    gene_col <- gene_col[!grepl("ID", gene_col)][1]
    out <- data.frame(probe = tab[[probe_col]], gene = trimws(tab[[gene_col]]),
                      stringsAsFactors = FALSE)
  }
  out <- out[!is.na(out$gene) & out$gene != "" & out$gene != "---", ]
  out <- out[!grepl("///|^\\?", out$gene), ]
  # Handle "A // B" multi-gene annotations: take the first
  out$gene <- sub(" //.*$", "", out$gene)
  out <- out[!duplicated(out$probe), ]
  out
}

# ---------- Cohort configuration ----------
cohorts <- list(
  GSE20347 = list(
    file = "GSE20347_series_matrix.txt.gz",
    anno = "GPL571.annot.gz", anno_type = "annot",
    normal_pattern = "^normal",
    log2 = TRUE
  ),
  GSE45670 = list(
    file = "GSE45670_series_matrix.txt.gz",
    anno = "GPL570.annot.gz", anno_type = "annot",
    normal_pattern = "^normal",
    log2 = FALSE   # linear values; log2 needed
  ),
  GSE75241 = list(
    file = "GSE75241_series_matrix.txt.gz",
    anno = "GPL5175_family.soft.gz", anno_type = "soft",
    normal_pattern = "^nonmalignant",
    log2 = TRUE
  )
  # GSE53624 is handled separately once downloaded (includes survival data)
)

# ---------- Main loop ----------
for (gse in names(cohorts)) {
  cfg <- cohorts[[gse]]
  gz <- file.path(geo_dir, cfg$file)
  ap <- file.path(anno_dir, cfg$anno)
  if (!file.exists(gz)) { cat(gse, ": matrix file missing; skipped\n"); next }
  if (!file.exists(ap)) { cat(gse, ": platform annotation missing; skipped\n"); next }
  cat("===== Processing", gse, "=====\n")
  sm <- read_series_matrix(gz)
  probe2gene <- parse_annot(ap, cfg$anno_type)
  cat("Platform probe-annotation rows:", nrow(probe2gene), "\n")
  expr <- sm$expr
  rownames(expr) <- expr[, 1]
  expr <- expr[, -1, drop = FALSE]
  # Convert columns to numeric (handling possible non-numeric characters; use [] <- to preserve data-frame structure and row names)
  expr[] <- lapply(expr, function(x) suppressWarnings(as.numeric(as.character(x))))
  # Keep annotated probes
  common <- intersect(rownames(expr), probe2gene$probe)
  cat("Probes:", nrow(expr), "-> annotated:", length(common), "\n")
  expr <- expr[common, , drop = FALSE]
  gene_expr <- collapse_probes(expr, probe2gene)
  cat("Genes after collapsing:", nrow(gene_expr), "\n")
  # Numeric transformation
  if (!cfg$log2) {
    gene_expr[gene_expr < 0] <- 0
    gene_expr <- log2(gene_expr + 1)
  }
  # Phenotype: tissue type
  titles <- sm$meta$title
  tissue <- ifelse(grepl(cfg$normal_pattern, titles, ignore.case = TRUE), "Normal", "Tumor")
  pheno <- data.frame(sample = colnames(gene_expr), gsm = sm$meta$geo_accession,
                      title = titles, tissue = tissue, stringsAsFactors = FALSE)
  cat("Tissue types:", paste(names(table(pheno$tissue)), table(pheno$tissue), sep = "=", collapse = ", "), "\n")
  saveRDS(list(expr = gene_expr, pheno = pheno,
               meta = list(gse = gse, platform = sm$platform, source = cfg$file)),
          file.path(out_dir, paste0("geo_", gse, ".rds")))
  cat("Saved geo_", gse, ".rds\n", sep = "")
}
cat("\n02_GEO_preprocess done\n")
