## ============================================================
## Figure organisation: results/*.png -> figures/ with unified naming + Fig2 merge
## ============================================================
suppressMessages({library(png); library(grid)})

rename <- function(from, to) {
  file.copy(from, to, overwrite = TRUE)
  cat(sprintf("%-40s -> %s\n", from, to))
}

rename("results/km_tcga.png",               "figures/fig6_km.png")
rename("results/roc_curves_tcga.png",       "figures/fig7_roc.png")
rename("results/nomogram.png",              "figures/fig8_nomogram.png")
rename("results/enrich_dotplot_GO.png",     "figures/fig9_go.png")
rename("results/enrich_barplot_KEGG.png",   "figures/fig9b_kegg.png")
rename("results/GSEA_plots.png",            "figures/fig10_gsea_kegg.png")
rename("results/GSEA_hallmark_plots.png",   "figures/fig10b_gsea_hallmark.png")
rename("results/immune_boxplot.png",        "figures/fig11_immune_box.png")
rename("results/immune_cor_risk.png",       "figures/fig11b_immune_cor.png")
rename("results/immune_heatmap.png",        "figures/fig11c_immune_heatmap.png")

# ---- Fig2 merge: 2x2 (three volcano plots + Venn) ----
panel <- function(path) {
  img <- readPNG(path)
  rasterGrob(img, interpolate = TRUE)
}
# Draw a grob inside a layout cell (equivalent to print(grob, vp = cell);
# the print(grob, vp = ...) form renders blank on grid >= R 4.6)
draw_in_cell <- function(g, row, col) {
  pushViewport(viewport(layout.pos.row = row, layout.pos.col = col))
  grid.draw(g)
  popViewport()
}

png("figures/fig2_deg_panel.png", width = 2600, height = 2400, res = 200)
pushViewport(viewport(layout = grid.layout(2, 2, widths = c(0.5, 0.5), heights = c(0.5, 0.5),
                                           respect = FALSE)))
draw_in_cell(panel("figures/fig2_volcano_GSE20347.png"), 1, 1)
draw_in_cell(panel("figures/fig2_volcano_GSE45670.png"), 1, 2)
draw_in_cell(panel("figures/fig2_volcano_TCGA.png"), 2, 1)
draw_in_cell(panel("figures/fig2_venn.png"), 2, 2)
popViewport()
dev.off()
cat("figures/fig2_deg_panel.png generated\n")

cat("\n=== figures directory listing ===\n")
for (f in list.files("figures", pattern = "png$", full.names = TRUE)) {
  cat(sprintf("%-40s %s bytes\n", basename(f), format(file.size(f), big.mark = ",")))
}
