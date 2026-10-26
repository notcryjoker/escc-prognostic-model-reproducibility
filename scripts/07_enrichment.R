# ============================================================================
# 07_enrichment.R
# Purpose: functional enrichment analysis (GO / KEGG) of consensus DEGs + genome-wide GSEA
# Input: results/DEG_consensus.rds, results/DEG_TCGA.tsv (for ranking)
# Output: results/enrich_GO_*.tsv, enrich_KEGG.tsv, GSEA_KEGG_*.tsv
#       results/enrich_dotplot_GO.png, enrich_barplot_KEGG.png,
#       results/GSEA_plots.png
# ============================================================================
rm(list = ls())
suppressPackageStartupMessages({
  library(clusterProfiler); library(org.Hs.eg.db); library(DOSE)
  library(data.table); library(ggplot2)
})

base <- Sys.getenv("ESCC_ROOT", unset = "")

if (!nzchar(base)) {

  base <- tryCatch(dirname(dirname(normalizePath(sys.frame(1)$ofile))),

                error = function(e) getwd())

}
out_dir <- file.path(base, "results")

consensus <- readRDS(file.path(out_dir, "DEG_consensus.rds"))
cons_degs <- consensus$gene[consensus$is_consensus]
cat("Consensus DEGs:", length(cons_degs), "\n")

# Gene symbols -> ENTREZ IDs (required for KEGG)
gene_entrez <- bitr(cons_degs, fromType = "SYMBOL", toType = "ENTREZID",
                    OrgDb = org.Hs.eg.db)
cat("Successfully mapped to ENTREZ:", nrow(gene_entrez), "/", length(cons_degs), "\n")
entrez_ids <- gene_entrez$ENTREZID

# ---------- 1. GO enrichment ----------
go_bp <- enrichGO(gene = entrez_ids, OrgDb = org.Hs.eg.db, ont = "BP",
                  pAdjustMethod = "BH", qvalueCutoff = 0.05,
                  readable = TRUE)
go_mf <- enrichGO(gene = entrez_ids, OrgDb = org.Hs.eg.db, ont = "MF",
                  pAdjustMethod = "BH", qvalueCutoff = 0.05,
                  readable = TRUE)
go_cc <- enrichGO(gene = entrez_ids, OrgDb = org.Hs.eg.db, ont = "CC",
                  pAdjustMethod = "BH", qvalueCutoff = 0.05,
                  readable = TRUE)
for (nm in c("BP", "MF", "CC")) {
  obj <- switch(nm, BP = go_bp, MF = go_mf, CC = go_cc)
  cat("GO-", nm, " significant terms:", nrow(obj), "\n", sep = "")
  if (nrow(obj) > 0)
    fwrite(as.data.frame(obj), file.path(out_dir, paste0("enrich_GO_", nm, ".tsv")), sep = "\t")
}

# ---------- 2. KEGG enrichment (local MSigDB GMT; KEGG REST unreachable) ----------
msigdb_dir <- file.path(base, "data/msigdb")
read_gmt <- function(f) {
  lines <- readLines(f)
  genes <- unlist(lapply(strsplit(lines, "\t"), function(x) x[-(1:2)]))
  terms <- unlist(lapply(strsplit(lines, "\t"), function(x) rep(x[1], length(x) - 2)))
  data.frame(term = terms, gene = genes, stringsAsFactors = FALSE)
}
t2g_kegg <- read_gmt(file.path(msigdb_dir, "c2.cp.kegg.v2023.1.Hs.entrez.gmt"))
cat("KEGG gene-set terms:", length(unique(t2g_kegg$term)), "\n")

kegg <- enricher(gene = entrez_ids, TERM2GENE = t2g_kegg,
                 pAdjustMethod = "BH", qvalueCutoff = 0.05)
kegg <- setReadable(kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
cat("KEGG significant terms:", nrow(kegg), "\n")
if (nrow(kegg) > 0)
  fwrite(as.data.frame(kegg), file.path(out_dir, "enrich_KEGG.tsv"), sep = "\t")

# ---------- 3. GSEA (all genes ranked by TCGA logFC) ----------
tcga_deg <- fread(file.path(out_dir, "DEG_TCGA.tsv"), data.table = FALSE)
tcga_deg <- tcga_deg[!is.na(tcga_deg$logFC), ]
tcga_deg <- tcga_deg[order(-tcga_deg$logFC), ]
# Symbol -> ENTREZ mapping (deduplicate: keep the largest |logFC|)
tcga_map <- bitr(unique(tcga_deg$gene), fromType = "SYMBOL", toType = "ENTREZID",
                 OrgDb = org.Hs.eg.db)
tcga_deg <- merge(tcga_deg, tcga_map, by.x = "gene", by.y = "SYMBOL")
tcga_deg <- tcga_deg[order(-abs(tcga_deg$logFC)), ]
tcga_deg <- tcga_deg[!duplicated(tcga_deg$ENTREZID), ]
gene_list <- sort(setNames(tcga_deg$logFC, tcga_deg$ENTREZID), decreasing = TRUE)
cat("Genes ranked for GSEA:", length(gene_list), "\n")

# KEGG GSEA (local gene sets)
gsea_kegg <- try(GSEA(geneList = gene_list, TERM2GENE = t2g_kegg,
                      pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,
                      pvalueCutoff = 0.05, verbose = FALSE, seed = TRUE), silent = TRUE)
if (!inherits(gsea_kegg, "try-error") && nrow(gsea_kegg) > 0) {
  gsea_kegg <- setReadable(gsea_kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  fwrite(as.data.frame(gsea_kegg), file.path(out_dir, "GSEA_KEGG.tsv"), sep = "\t")
  cat("GSEA-KEGG significant terms:", nrow(gsea_kegg), "\n")
} else {
  cat("GSEA-KEGG: no significant terms or failed\n")
}

# Hallmark GSEA
t2g_hall <- read_gmt(file.path(msigdb_dir, "h.all.v2024.1.Hs.entrez.gmt"))
gsea_hall <- try(GSEA(geneList = gene_list, TERM2GENE = t2g_hall,
                      pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,
                      pvalueCutoff = 0.05, verbose = FALSE, seed = TRUE), silent = TRUE)
if (!inherits(gsea_hall, "try-error") && nrow(gsea_hall) > 0) {
  gsea_hall <- setReadable(gsea_hall, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  fwrite(as.data.frame(gsea_hall), file.path(out_dir, "GSEA_hallmark.tsv"), sep = "\t")
  cat("GSEA-Hallmark significant terms:", nrow(gsea_hall), "\n")
} else {
  cat("GSEA-Hallmark: no significant terms or failed\n")
}

# ---------- 4. Visualisation (custom ggplot2, avoids version-compatibility issues) ----------
# 4.1 GO-BP dot plot (top 10)
if (exists("go_bp") && !is.null(go_bp) && nrow(go_bp) > 0) {
  d <- as.data.frame(go_bp)
  d <- d[order(d$p.adjust), ]
  d <- d[1:min(10, nrow(d)), ]
  d$GeneRatio <- sapply(strsplit(d$GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
  p <- ggplot(d, aes(x = GeneRatio, y = reorder(Description, GeneRatio),
                     color = p.adjust, size = Count)) +
    geom_point(alpha = 0.85) +
    scale_color_gradient(low = "#B2182B", high = "#2166AC",
                         name = "Adj. P", trans = "log10") +
    scale_size_continuous(range = c(3, 8), name = "Gene count") +
    labs(x = "Gene ratio", y = NULL,
         title = "GO Biological Process enrichment (top 10)") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.y = element_text(size = 9))
  ggsave(file.path(out_dir, "enrich_dotplot_GO.png"), p, width = 9, height = 7, dpi = 300)
}

# 4.2 KEGG bar plot (top 15)
if (exists("kegg") && !is.null(kegg) && nrow(kegg) > 0) {
  d <- as.data.frame(kegg)
  d <- d[order(d$p.adjust), ]
  d <- d[1:min(15, nrow(d)), ]
  p <- ggplot(d, aes(x = Count, y = reorder(Description, Count),
                     fill = p.adjust)) +
    geom_col(alpha = 0.9) +
    scale_fill_gradient(low = "#B2182B", high = "#2166AC",
                        name = "Adj. P", trans = "log10") +
    labs(x = "Gene count", y = NULL, title = "KEGG pathway enrichment (top 15)") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.y = element_text(size = 9))
  ggsave(file.path(out_dir, "enrich_barplot_KEGG.png"), p, width = 9, height = 7, dpi = 300)
}

# 4.3 GSEA running-score curves (custom, top 2 significant KEGG pathways)
plot_running_score <- function(gene_list, geneset, title) {
  gl <- gene_list[!is.na(gene_list)]
  hit <- names(gl) %in% geneset
  n <- length(gl); nh <- sum(hit)
  inc <- ifelse(hit, abs(gl), 0); dec <- ifelse(hit, 0, 1)
  # Normalisation: hits weighted by |stat|, misses equally weighted
  sum_hit <- sum(inc)
  if (sum_hit == 0) return(NULL)
  run <- cumsum(inc / sum_hit - dec / (n - nh))
  run <- run - min(run)   # start from 0
  df <- data.frame(pos = seq_along(run), score = run, hit = hit)
  max_i <- which.max(df$score); min_i <- which.min(df$score)
  es <- if (abs(df$score[max_i]) >= abs(df$score[min_i])) df$score[max_i] else df$score[min_i]
  p <- ggplot(df, aes(x = pos, y = score)) +
    geom_line(color = "#2166AC", linewidth = 0.8) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
    geom_vline(xintercept = df$pos[df$hit], color = "#B2182B",
               alpha = 0.35, linewidth = 0.5) +
    annotate("text", x = max(df$pos) * 0.98, y = max(df$score) * 0.92,
             label = sprintf("ES = %.3f", es), hjust = 1, size = 3.5, fontface = "bold") +
    labs(x = "Gene rank (ordered by log2FC)", y = "Running enrichment score",
         title = title) +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 11))
  p
}

if (exists("gsea_kegg") && !inherits(gsea_kegg, "try-error") && nrow(gsea_kegg) >= 1) {
  top_sets <- gsea_kegg$ID[1:min(2, nrow(gsea_kegg))]
  plots <- lapply(top_sets, function(sid) {
    gs <- t2g_kegg$gene[t2g_kegg$term == sid]
    plot_running_score(gene_list, gs, sid)
  })
  plots <- plots[!sapply(plots, is.null)]
  if (length(plots) > 0)
    ggsave(file.path(out_dir, "GSEA_plots.png"),
           patchwork::wrap_plots(plots, ncol = 1) |> (\(p) p)(),
           width = 9, height = 4.5 * length(plots), dpi = 300)
}

if (exists("gsea_hall") && !inherits(gsea_hall, "try-error") && nrow(gsea_hall) >= 1) {
  top_sets <- gsea_hall$ID[1:min(2, nrow(gsea_hall))]
  plots <- lapply(top_sets, function(sid) {
    gs <- t2g_hall$gene[t2g_hall$term == sid]
    plot_running_score(gene_list, gs, sid)
  })
  plots <- plots[!sapply(plots, is.null)]
  if (length(plots) > 0)
    ggsave(file.path(out_dir, "GSEA_hallmark_plots.png"),
           patchwork::wrap_plots(plots, ncol = 1) |> (\(p) p)(),
           width = 9, height = 4.5 * length(plots), dpi = 300)
}

cat("\n07_enrichment done\n")
