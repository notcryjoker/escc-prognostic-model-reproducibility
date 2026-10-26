# ============================================================================
# 01b_TCGA_to_symbol.R
# Purpose: convert the TCGA expression matrix from ENSG IDs to gene symbols (harmonised with the GEO cohorts)
# Method: use the gene_name column of the STAR-counts files; for one symbol with several ENSGs keep the max TPM
# Input: data/tcga/star_counts/*.tsv, results/tcga_expr_tpm.rds
# Output: results/tcga_expr_symbol.rds    tumour-sample symbol matrix (log2(TPM+1))
#        results/tcga_normal_symbol.rds normal-sample symbol matrix (log2(TPM+1))
# ============================================================================
rm(list = ls())
suppressPackageStartupMessages({library(data.table); library(matrixStats)})

base <- Sys.getenv("ESCC_ROOT", unset = "")

if (!nzchar(base)) {

  base <- tryCatch(dirname(dirname(normalizePath(sys.frame(1)$ofile))),

                error = function(e) getwd())

}
out_dir <- file.path(base, "results")

# ---------- 1. Build the ENSG -> symbol mapping ----------
fp <- list.files(file.path(base, "data/tcga/star_counts"), full.names = TRUE)[1]
anno <- fread(fp, data.table = FALSE, showProgress = FALSE)
anno <- anno[!grepl("^N_", anno$gene_id), ]
anno$gene_id <- sub("\\.[0-9]+$", "", anno$gene_id)
ensg2sym <- setNames(anno$gene_name, anno$gene_id)
cat("Mapping entries:", length(ensg2sym), "| unique symbols:", length(unique(ensg2sym)), "\n")

to_symbol <- function(mat, label) {
  rn <- sub("\\.[0-9]+$", "", rownames(mat))
  sym <- ensg2sym[rn]
  keep <- !is.na(sym) & sym != "" & !grepl("^NA$", sym)
  cat(label, "mappable genes:", sum(keep), "/", nrow(mat), "\n")
  m <- mat[keep, ]
  s <- sym[keep]
  grp <- split(seq_len(nrow(m)), s)
  out <- do.call(rbind, lapply(grp, function(idx) {
    if (length(idx) == 1) m[idx, ] else colMaxs(m[idx, , drop = FALSE])
  }))
  rownames(out) <- names(grp)
  cat(label, "symbol matrix dimensions:", dim(out), "\n")
  out
}

# ---------- 2. Tumour matrix ----------
tc <- readRDS(file.path(out_dir, "tcga_expr_tpm.rds"))
expr_sym <- to_symbol(tc$expr, "[tumour]")
saveRDS(expr_sym, file.path(out_dir, "tcga_expr_symbol.rds"))

# ---------- 3. Normal matrix (3 Solid Tissue Normal samples) ----------
fj <- jsonlite::fromJSON(file.path(base, "data/tcga/ESCC_rnaseq_filelist.json"),
                         simplifyVector = FALSE)
normal_files <- sapply(fj, function(h)
  if (h$cases[[1]]$samples[[1]]$sample_type == "Solid Tissue Normal") h$file_id else NA_character_)
normal_files <- normal_files[!is.na(normal_files)]
cat("Normal-sample files:", length(normal_files), "\n")

normal_expr <- list()
for (fid in normal_files) {
  fp <- file.path(base, "data/tcga/star_counts", paste0(fid, ".tsv"))
  if (!file.exists(fp)) { cat("Missing:", fid, "\n"); next }
  d <- fread(fp, data.table = FALSE, showProgress = FALSE)
  d <- d[!grepl("^N_", d$gene_id), ]
  d$gene_id <- sub("\\.[0-9]+$", "", d$gene_id)
  case <- sapply(fj, function(h) if (h$file_id == fid) h$cases[[1]]$submitter_id else NULL)
  case <- unlist(case)[1]
  normal_expr[[case]] <- setNames(log2(d$tpm_unstranded + 1), d$gene_id)
}
ne <- do.call(cbind, normal_expr)
normal_sym <- to_symbol(ne, "[normal]")
saveRDS(normal_sym, file.path(out_dir, "tcga_normal_symbol.rds"))

# ---------- 4. Overlap check with the GEO cohorts ----------
g1 <- readRDS(file.path(out_dir, "geo_GSE20347.rds"))
g2 <- readRDS(file.path(out_dir, "geo_GSE45670.rds"))
cat("Overlap with GSE20347:", length(intersect(rownames(expr_sym), rownames(g1$expr))), "\n")
cat("Overlap with GSE45670:", length(intersect(rownames(expr_sym), rownames(g2$expr))), "\n")
cat("\n01b_TCGA_to_symbol done\n")
