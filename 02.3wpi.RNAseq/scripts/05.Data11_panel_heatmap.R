#!/usr/bin/env Rscript
# ============================================================================
# 05.Data11_panel_heatmap.R
# Reproduce the J002 "Data 11" curated nodulation-gene heatmap on J003 data.
# ----------------------------------------------------------------------------
# The gene panel (41 genes, with their display labels and functional-block
# order) was transcribed from the original figure and lives in:
#     Data11_panel_genes.tsv   (columns: category, label, suffix)
# where `suffix` is the locus tail (e.g. Chr5g0427351) appended to "MtrunA17"
# to form the gene_name used in MasterTable_with_DE.tsv. All 41 verified
# present in the J003 filtered set.
#
# GREY (NA) SEMANTICS -- this is the key reproduction decision:
#   The original masks log2FC by significance (the entire rsp86 column is grey
#   because almost nothing is a DEG there). We reproduce that as the primary
#   figure (MASK_BY_SIG = significance-masked), and ALSO emit an all-log2FC
#   variant where every tested cell is colored.
#
# COLUMNS: the 5 A-genotypes in the original order (A86, A256, A485, A86256,
#   A86485) plus dnf2 as a 6th reference column.
#
# Row order and functional blocks are FIXED (no clustering) to match the
# original. Colors: blue (down) - white - red (up), centered at 0, clipped +/-5.
#
# USAGE: Rscript 05.Data11_panel_heatmap.R
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(pheatmap)
  library(RColorBrewer)
  library(grid)
})

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
BASE_DIR <- ".."
DE_DIR   <- file.path(BASE_DIR, "05.DE_analysis")
SUM_DIR  <- file.path(DE_DIR, "99.Summary")
MASTER   <- file.path(SUM_DIR, "MasterTable_with_DE.tsv")
PANEL    <- file.path(DE_DIR, "Data11_panel_genes.tsv")

# Columns in the original's order, + dnf2.
GENOTYPES   <- c("A86", "A256", "A485", "A86256", "A86485", "dnf2")
# Pretty column headers echoing the original rsp* naming.
COL_LABELS  <- c("rsp86", "rsp256", "rsp485", "rsp86/256", "rsp86/485", "dnf2")

LFC_CLIP    <- 5          # symmetric color clip (data range is [-4.65, 5.65])
IMG_FORMATS <- c("pdf", "png", "tiff")
IMG_DPI     <- 600

# ---------------------------------------------------------------------------
# Save a pheatmap object to all formats (silent=TRUE returns a gtable that
# must be grid.draw'n explicitly, else blank pages).
# ---------------------------------------------------------------------------
save_pheatmap <- function(ph, base, width, height) {
  for (fmt in IMG_FORMATS) {
    fp <- file.path(SUM_DIR, paste0(base, ".", fmt))
    if (fmt == "pdf")  pdf(fp, width = width, height = height)
    if (fmt == "png")  png(fp, width = width, height = height, units = "in", res = IMG_DPI)
    if (fmt == "tiff") tiff(fp, width = width, height = height, units = "in",
                            res = IMG_DPI, compression = "lzw")
    grid::grid.newpage(); grid::grid.draw(ph$gtable); dev.off()
  }
}

# ---------------------------------------------------------------------------
# Load data + panel; resolve every gene by gene_name (NOT gene_id).
# ---------------------------------------------------------------------------
cat("== Loading master table and panel ==\n")
m     <- fread(MASTER)
panel <- fread(PANEL)
panel[, gene_name := paste0("MtrunA17", suffix)]

idx <- match(panel$gene_name, m$gene_name)
if (anyNA(idx))
  stop("panel genes not found in master: ",
       paste(panel$label[is.na(idx)], collapse = ", "))
cat(sprintf("   resolved %d / %d panel genes\n", sum(!is.na(idx)), nrow(panel)))

lfc_cols <- paste0(GENOTYPES, "_logFC")
sig_cols <- paste0(GENOTYPES, "_sig")
fdr_cols <- paste0(GENOTYPES, "_FDR")

mrows <- m[idx]                              # 41 rows, panel order preserved
stopifnot(identical(mrows$gene_name, panel$gene_name))   # integrity: same order

# log2FC matrix in panel row-order, genotype col-order.
lfc <- as.matrix(mrows[, ..lfc_cols])
rownames(lfc) <- panel$label
colnames(lfc) <- COL_LABELS

# Significance mask (TRUE where that genotype called the gene a DEG).
sig <- as.matrix(mrows[, ..sig_cols]) != "NS"
dimnames(sig) <- dimnames(lfc)

# Significance stars from FDR (only where significant).
fdr <- as.matrix(mrows[, ..fdr_cols])
star <- matrix("", nrow(lfc), ncol(lfc), dimnames = dimnames(lfc))
star[sig] <- ifelse(fdr[sig] < 1e-3, "***", ifelse(fdr[sig] < 1e-2, "**", "*"))

# ---------------------------------------------------------------------------
# Row annotation: functional category blocks + gap rows between blocks.
# Category order = first appearance in the panel file (preserves the figure).
# ---------------------------------------------------------------------------
cat_levels <- unique(panel$category)
row_anno   <- data.frame(Category = factor(panel$category, levels = cat_levels),
                         row.names = panel$label)
cat_pal    <- setNames(colorRampPalette(brewer.pal(8, "Set2"))(length(cat_levels)),
                       cat_levels)
gaps_rows  <- head(cumsum(table(factor(panel$category, levels = cat_levels))), -1)

# Diverging palette, centered at 0, clipped to +/- LFC_CLIP.
breaks  <- seq(-LFC_CLIP, LFC_CLIP, length.out = 101)
pal     <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100)
clip    <- function(x) pmax(pmin(x, LFC_CLIP), -LFC_CLIP)

draw_panel <- function(mat, base, title) {
  ph <- pheatmap(
    clip(mat),
    color = pal, breaks = breaks, scale = "none",
    cluster_rows = FALSE, cluster_cols = FALSE,
    gaps_row = gaps_rows,
    annotation_row = row_anno, annotation_colors = list(Category = cat_pal),
    annotation_names_row = FALSE,
    display_numbers = star, number_color = "black", fontsize_number = 6,
    cellwidth = 24, cellheight = 11,
    fontsize_row = 7, fontsize_col = 9, angle_col = 45,
    border_color = "grey85", na_col = "grey85",
    legend_breaks = c(-5, -2.5, 0, 2.5, 5),
    legend_labels = c("-5", "-2.5", "0", "+2.5", "+5"),
    main = title, silent = TRUE)
  save_pheatmap(ph, base, width = 7.5, height = 11)
}

# ---------------------------------------------------------------------------
# (1) Significance-masked = faithful reproduction of the original.
# (2) All-log2FC = magnitude everywhere a value exists.
# ---------------------------------------------------------------------------
cat("== Drawing significance-masked panel (reproduction) ==\n")
lfc_masked <- lfc; lfc_masked[!sig] <- NA
draw_panel(lfc_masked, "Data11_panel_heatmap_sigmasked",
           "Curated symbiosis panel - log2FC vs WT (significant only)")

cat("== Drawing all-log2FC panel ==\n")
draw_panel(lfc, "Data11_panel_heatmap_alllogFC",
           "Curated symbiosis panel - log2FC vs WT (all values)")

# ---------------------------------------------------------------------------
# Export the underlying data table (log2FC + FDR + sig for all 6 genotypes).
# ---------------------------------------------------------------------------
out <- cbind(
  data.table(label = panel$label, category = panel$category,
             gene_name = panel$gene_name,
             gene_id = mrows$gene_id, acronym = mrows$acronym,
             geneProduct = mrows$geneProduct),
  as.data.table(setNames(as.data.frame(lfc), paste0(GENOTYPES, "_logFC"))),
  as.data.table(setNames(as.data.frame(fdr), paste0(GENOTYPES, "_FDR"))),
  as.data.table(setNames(as.data.frame(mrows[, ..sig_cols]), paste0(GENOTYPES, "_sig")))
)
fwrite(out, file.path(SUM_DIR, "Data11_panel_heatmap_data.tsv"), sep = "\t")

cat("\n== DONE ==\n")
cat("Genes:", nrow(panel), "| Genotypes:", paste(GENOTYPES, collapse=", "), "\n")
cat("Outputs under:", SUM_DIR, "\n")
