## ============================================================
## Export supplementary tables as BMC "Additional file N" (see manuscript supplementary section)
## ============================================================
suppressMessages({ library(data.table) })
out_dir <- "results"
sup_dir <- "manuscript/supplementary"
dir.create(sup_dir, showWarnings = FALSE, recursive = TRUE)

# ---- S1: 719 consensus DEGs (exported with the original column names) ----
deg <- fread(file.path(out_dir, "DEG_consensus.tsv"))
fwrite(deg, file.path(sup_dir, "Additional file 4.tsv"), sep = "\t")
cat("S1 consensus DEGs:", nrow(deg), "rows, columns:", paste(names(deg), collapse = ", "), "\n")

# ---- S2: 48 univariable Cox candidate genes ----
pg <- fread(file.path(out_dir, "prognostic_genes.tsv"))
fwrite(pg, file.path(sup_dir, "Additional file 5.tsv"), sep = "\t")
cat("S2 candidate genes:", nrow(pg), "rows -> TableS2_univariate_cox_candidates.tsv\n")

# ---- S3a: GO enrichment ----
go <- rbind(
  cbind(Ontology = "BP", fread(file.path(out_dir, "enrich_GO_BP.tsv"))),
  cbind(Ontology = "CC", fread(file.path(out_dir, "enrich_GO_CC.tsv"))),
  cbind(Ontology = "MF", fread(file.path(out_dir, "enrich_GO_MF.tsv"))),
  fill = TRUE)
fwrite(go, file.path(sup_dir, "Additional file 7.tsv"), sep = "\t")
cat("S3a GO:", nrow(go), "rows\n")

# ---- S3b: KEGG ORA ----
kegg <- fread(file.path(out_dir, "enrich_KEGG.tsv"))
fwrite(kegg, file.path(sup_dir, "Additional file 8.tsv"), sep = "\t")
cat("S3b KEGG ORA:", nrow(kegg), "rows\n")

# ---- S3c/d: GSEA ----
gk <- fread(file.path(out_dir, "GSEA_KEGG.tsv"))
gh <- fread(file.path(out_dir, "GSEA_hallmark.tsv"))
fwrite(gk, file.path(sup_dir, "Additional file 9.tsv"), sep = "\t")
fwrite(gh, file.path(sup_dir, "Additional file 10.tsv"), sep = "\t")
cat("S3c GSEA-KEGG:", nrow(gk), "rows | S3d GSEA-Hallmark:", nrow(gh), "rows\n")

# ---- S4: immune infiltration ----
imm <- fread(file.path(out_dir, "immune_high_vs_low.tsv"))
cor <- fread(file.path(out_dir, "immune_cor_risk.tsv"))
fwrite(imm, file.path(sup_dir, "Additional file 11.tsv"), sep = "\t")
fwrite(cor, file.path(sup_dir, "Additional file 12.tsv"), sep = "\t")
cat("S4 immune:", nrow(imm), "cell types\n")

# ---- S5: immune-checkpoint correlations ----
ck <- fread(file.path(out_dir, "immune_checkpoint_cor.tsv"))
fwrite(ck, file.path(sup_dir, "Additional file 13.tsv"), sep = "\t")
cat("S5 immune checkpoints:", nrow(ck), "\n")

# ---- S6: split validation ----
sv <- fread(file.path(out_dir, "split_validation.tsv"))
svf <- fread(file.path(out_dir, "split_validation_fixedmodel.tsv"))
fwrite(sv, file.path(sup_dir, "Additional file 14.tsv"), sep = "\t")
fwrite(svf, file.path(sup_dir, "Additional file 15.tsv"), sep = "\t")
cat("S6 split validation:", nrow(sv), "/", nrow(svf), "runs\n")

# ---- S7: 17-gene model coefficients ----
mdl <- readRDS(file.path(out_dir, "lasso_model.rds"))
coef_df <- data.frame(gene = mdl$genes,
                      coefficient = round(as.numeric(mdl$coef), 4),
                      role = ifelse(as.numeric(mdl$coef) > 0, "risk", "protective"))
coef_df <- coef_df[order(-coef_df$coefficient), ]
fwrite(coef_df, file.path(sup_dir, "Additional file 6.tsv"), sep = "\t")
cat("S7 model coefficients:", nrow(coef_df), "genes\n")

cat("\nAll supplementary materials exported to:", normalizePath(sup_dir), "\n")
