# ============================================================================
# 01_TCGA_preprocess.R
# Purpose: merge TCGA-ESCC STAR-counts expression data; build the expression matrix and clinical/survival table
# Input: data/tcga/star_counts/*.tsv (GDC augmented_star_gene_counts)
#       data/tcga/ESCC_clinical.tsv
# Output: results/tcga_expr_tpm.rds      log2(TPM+1) expression matrix
#        results/tcga_counts.rds        unstranded raw count matrix
#        results/tcga_clinical_surv.rds clinical + survival information (matched to analysis samples)
#        results/tcga_qc_stats.txt     QC statistics
# ============================================================================
rm(list = ls())
suppressPackageStartupMessages({
  library(data.table); library(matrixStats)
})

base <- Sys.getenv("ESCC_ROOT", unset = "")

if (!nzchar(base)) {

  base <- tryCatch(dirname(dirname(normalizePath(sys.frame(1)$ofile))),

                error = function(e) getwd())

}
star_dir <- file.path(base, "data/tcga/star_counts")
out_dir  <- file.path(base, "results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- 1. Read the file manifest (file_id -> case_id) ----
fl <- fread(file.path(base, "data/tcga/download_list.txt"), header = FALSE,
            col.names = c("file_id", "case_id"))
# Take sample_type from the file-manifest JSON; keep Primary Tumor only
fj <- jsonlite::fromJSON(file.path(base, "data/tcga/ESCC_rnaseq_filelist.json"),
                         simplifyVector = FALSE)
sample_info <- data.frame(
  file_id = sapply(fj, function(h) h$file_id),
  case_id = sapply(fj, function(h) h$cases[[1]]$submitter_id),
  sample_id = sapply(fj, function(h) h$cases[[1]]$samples[[1]]$submitter_id),
  sample_type = sapply(fj, function(h) h$cases[[1]]$samples[[1]]$sample_type),
  stringsAsFactors = FALSE
)
cat("File manifest:", nrow(sample_info), "files\n")
print(table(sample_info$sample_type))

# ---- 2. Read files one by one; extract unstranded counts and TPM ----
files <- list.files(star_dir, pattern = "\\.tsv$", full.names = TRUE)
cat("Files on disk:", length(files), "\n")

count_list <- list(); tpm_list <- list(); meta_list <- list()
for (fp in files) {
  fid <- sub("\\.tsv$", "", basename(fp))
  si <- sample_info[sample_info$file_id == fid, ]
  if (nrow(si) == 0) next
  # Keep Primary Tumor only
  si <- si[si$sample_type == "Primary Tumor", ]
  if (nrow(si) == 0) next
  d <- fread(fp, header = TRUE, data.table = FALSE, showProgress = FALSE)
  # The first 4-5 rows of GDC augmented STAR files are summary statistics (N_unmapped, etc.)
  d <- d[!grepl("^N_", d$gene_id), ]
  colnames(d)[1:3] <- c("gene_id", "gene_name", "gene_type")
  # Remove ENSEMBL gene ID version suffixes (ENSG00000000003.15 -> ENSG00000000003)
  d$gene_id <- sub("\\.[0-9]+$", "", d$gene_id)
  gene_key <- d$gene_id
  count_list[[si$case_id[1]]] <- setNames(d$unstranded, gene_key)
  tpm_list[[si$case_id[1]]]   <- setNames(d$tpm_unstranded, gene_key)
  meta_list[[si$case_id[1]]]  <- si[1, ]
}

cat("Samples successfully read:", length(count_list), "\n")

# If a case has multiple files (keep the one with the higher median count)
dup_cases <- names(which(table(names(count_list)) > 1))
if (length(dup_cases) > 0) {
  cat("Note: duplicate cases:", paste(dup_cases, collapse = ", "), "\n")
}

# ---- 3. Build matrices ----
all_genes <- unique(unlist(lapply(count_list, names)))
count_mat <- do.call(cbind, lapply(count_list, function(x) x[all_genes]))
tpm_mat   <- do.call(cbind, lapply(tpm_list, function(x) x[all_genes]))
colnames(count_mat) <- colnames(tpm_mat) <- names(count_list)
rownames(count_mat) <- rownames(tpm_mat) <- all_genes
# Gene annotation (taken from the first file)
d0 <- fread(files[1], header = TRUE, data.table = FALSE)
d0 <- d0[!grepl("^N_", d0$gene_id), ]
d0$gene_id <- sub("\\.[0-9]+$", "", d0$gene_id)
gene_anno <- d0[, c("gene_id", "gene_name", "gene_type")]
gene_anno <- gene_anno[match(all_genes, gene_anno$gene_id), ]
rownames(gene_anno) <- NULL

cat("Count matrix dimensions:", dim(count_mat), "\n")

# ---- 4. Filter low-expression genes ----
# Keep: TPM >= 1 or count >= 6 in >=50% of tumour samples
keep <- rowSums(tpm_mat >= 1) >= ceiling(0.5 * ncol(tpm_mat)) |
        rowSums(count_mat >= 6) >= ceiling(0.5 * ncol(count_mat))
cat("Genes retained after expression filtering:", sum(keep), "/", nrow(count_mat), "\n")
count_mat_f <- count_mat[keep, ]
tpm_mat_f   <- tpm_mat[keep, ]
gene_anno_f <- gene_anno[keep, ]

# ---- 5. Merge clinical and survival data ----
clin <- fread(file.path(base, "data/tcga/ESCC_clinical.tsv"), data.table = FALSE)
clin <- clin[clin$case_id %in% colnames(tpm_mat_f), ]
clin <- clin[match(colnames(tpm_mat_f), clin$case_id), ]
rownames(clin) <- clin$case_id

# Survival object: time = time to death or last follow-up; status = 1 for death
clin$os_time <- ifelse(clin$vital_status == "Dead",
                       as.numeric(clin$days_to_death),
                       as.numeric(clin$days_to_last_follow_up))
clin$os_status <- ifelse(clin$vital_status == "Dead", 1, 0)
clin$os_time[clin$os_time <= 0] <- NA
clin$stage <- clin$ajcc_pathologic_stage
clin$stage_group <- ifelse(grepl("^Stage IV", clin$stage), "IV",
                   ifelse(grepl("^Stage III", clin$stage), "III",
                   ifelse(grepl("^Stage II", clin$stage), "II",
                   ifelse(grepl("^Stage I", clin$stage), "I", NA))))
clin$age <- as.numeric(clin$age_at_index)

# ---- 6. Expression-matrix normalisation ----
# log2(TPM+1), used for all downstream analyses and modelling
expr_log <- log2(tpm_mat_f + 1)

# ---- 7. QC statistics ----
qc <- c(
  paste0("Samples (Primary Tumor): ", ncol(expr_log)),
  paste0("Genes (after filtering): ", nrow(expr_log)),
  paste0("Death events: ", sum(clin$os_status, na.rm = TRUE)),
  paste0("Complete survival data (time+status): ", sum(complete.cases(clin$os_time, clin$os_status))),
  paste0("Stage available: ", sum(!is.na(clin$stage_group)), "/", nrow(clin)),
  paste0("Stage distribution: ", paste(names(table(clin$stage_group)), table(clin$stage_group), sep="=", collapse=", ")),
  paste0("Median follow-up time (days): ", median(clin$os_time, na.rm = TRUE)),
  paste0("Follow-up time range (days): ", paste(range(clin$os_time, na.rm = TRUE), collapse = "-"))
)
qc_text <- paste(qc, collapse = "\n")
writeLines(qc_text, file.path(out_dir, "tcga_qc_stats.txt"))
cat("\n=== TCGA QC ===\n", qc_text, "\n")

# ---- 8. Save ----
saveRDS(list(expr = expr_log, counts = count_mat_f, anno = gene_anno_f),
        file.path(out_dir, "tcga_expr_tpm.rds"))
saveRDS(clin, file.path(out_dir, "tcga_clinical_surv.rds"))
cat("\nSaved: tcga_expr_tpm.rds, tcga_clinical_surv.rds\n")
cat("01_TCGA_preprocess done\n")
