#!/usr/bin/env Rscript
# ============================================================================
# 02.summary_upset_heatmaps.R
# Cross-genotype summary visualizations, built from MasterTable_with_DE.tsv.
# ----------------------------------------------------------------------------
# Produces, under 99.Summary/:
#   1. UpSet plots (UP and DOWN) over ALL 6 genotypes vs WT -- the modern
#      replacement for >3-way Venns; shows every intersection without the
#      visual impossibility of a 6-circle Venn.
#   2. Data-driven summary heatmap: logFC of the union of top DEGs across
#      genotypes (no biology assumed), clustered.
#   3. Curated symbiosis-panel heatmap: hand-picked nodulation/senescence
#      genes (ported from the J002 figure), fixed row order by category,
#      with FDR significance stars.
#
# All overlap logic uses the per-genotype *_sig columns (UP / DOWN / NS) from
# the master table, so it is guaranteed consistent with the per-comparison
# DE_results tables (same filtered gene universe, same significance calls).
#
# USAGE: Rscript 02.summary_upset_heatmaps.R
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(UpSetR)
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

GENOTYPES   <- c("A86", "A256", "A485", "A86256", "A86485", "dnf2")
IMG_FORMATS <- c("pdf", "png", "tiff")
IMG_DPI     <- 600
TOP_N       <- 50      # top DEGs per genotype for the data-driven heatmap
LFC_CLIP    <- 4       # clip logFC for color scale in both heatmaps

dir.create(SUM_DIR, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# Device + saver helpers (mirrors 00.edgeR_unified.R)
# ---------------------------------------------------------------------------
open_dev <- function(fp, fmt, width, height) {
  if (fmt == "pdf")  pdf(fp, width = width, height = height)
  if (fmt == "png")  png(fp, width = width, height = height, units = "in", res = IMG_DPI)
  if (fmt == "tiff") tiff(fp, width = width, height = height, units = "in",
                          res = IMG_DPI, compression = "lzw")
}
# UpSetR draws to the active device on print(); wrap accordingly.
save_upset <- function(obj, base, width = 9, height = 6) {
  for (fmt in IMG_FORMATS) {
    open_dev(file.path(SUM_DIR, paste0(base, ".", fmt)), fmt, width, height)
    print(obj); dev.off()
  }
}
# pheatmap(silent=TRUE) returns a gtable that must be grid.draw'n explicitly.
save_pheatmap <- function(ph, base, width = 8, height = 10) {
  for (fmt in IMG_FORMATS) {
    open_dev(file.path(SUM_DIR, paste0(base, ".", fmt)), fmt, width, height)
    grid::grid.newpage(); grid::grid.draw(ph$gtable); dev.off()
  }
}

# ---------------------------------------------------------------------------
# Load master DE table
# ---------------------------------------------------------------------------
cat("== Loading master DE table ==\n")
m <- fread(MASTER)
stopifnot(all(paste0(GENOTYPES, "_sig")   %in% names(m)),
          all(paste0(GENOTYPES, "_logFC") %in% names(m)),
          all(paste0(GENOTYPES, "_FDR")   %in% names(m)))
cat(sprintf("   %d genes x %d genotypes\n", nrow(m), length(GENOTYPES)))

# Build UP / DOWN gene-id sets per genotype from the *_sig columns.
up_sets   <- lapply(GENOTYPES, function(g) m$gene_id[m[[paste0(g, "_sig")]] == "UP"])
down_sets <- lapply(GENOTYPES, function(g) m$gene_id[m[[paste0(g, "_sig")]] == "DOWN"])
names(up_sets) <- names(down_sets) <- GENOTYPES

# ===========================================================================
# 1. UpSet plots over all 6 genotypes (UP and DOWN separately)
# ===========================================================================
cat("== 1. UpSet plots ==\n")
make_upset <- function(sets, base, title_dir) {
  sizes <- sapply(sets, length)
  keep  <- names(sizes)[sizes > 0]          # fromList drops empty sets anyway
  if (length(keep) < 2) { cat("   <2 non-empty sets for", title_dir, "- skipping\n"); return() }
  obj <- upset(fromList(sets[keep]),
               sets = rev(keep), keep.order = TRUE,
               order.by = "freq", nintersects = 30,
               mainbar.y.label = sprintf("%s DEGs in intersection", title_dir),
               sets.x.label = sprintf("%s DEGs per genotype", title_dir),
               text.scale = c(1.4, 1.3, 1.2, 1.1, 1.3, 1.1),
               point.size = 2.6, line.size = 0.9,
               main.bar.color = ifelse(title_dir == "UP", "#B2182B", "#2166AC"),
               sets.bar.color = ifelse(title_dir == "UP", "#B2182B", "#2166AC"))
  save_upset(obj, base, width = 10, height = 6)
  cat(sprintf("   %s UpSet saved (%d non-empty sets)\n", title_dir, length(keep)))
}
make_upset(up_sets,   "UpSet_UP",   "UP")
make_upset(down_sets, "UpSet_DOWN", "DOWN")

# Also dump the full intersection breakdown as a table (every non-empty combo).
intersection_table <- function(sets, direction) {
  universe <- unique(unlist(sets))
  if (!length(universe)) return(NULL)
  memb <- sapply(sets, function(s) universe %in% s)   # genes x genotypes logical
  combo <- apply(memb, 1, function(r) paste(names(sets)[r], collapse = "&"))
  dt <- as.data.table(table(combo))[order(-N)]
  setnames(dt, c("Intersection", "n_genes"))
  dt[, Direction := direction]
  dt[]
}
ix <- rbindlist(list(intersection_table(up_sets, "UP"),
                     intersection_table(down_sets, "DOWN")))
fwrite(ix, file.path(SUM_DIR, "Intersection_counts.tsv"), sep = "\t")
cat("   intersection breakdown -> Intersection_counts.tsv\n")

# ===========================================================================
# 2. Data-driven summary heatmap (logFC of union of top DEGs)
# ===========================================================================
cat("== 2. Data-driven top-DEG heatmap ==\n")
lfc_cols <- paste0(GENOTYPES, "_logFC")
fdr_cols <- paste0(GENOTYPES, "_FDR")
sig_cols <- paste0(GENOTYPES, "_sig")

# For each genotype, rank its significant genes by |logFC| and take top N.
top_ids <- unique(unlist(lapply(GENOTYPES, function(g) {
  sub <- m[get(paste0(g, "_sig")) != "NS"]
  sub <- sub[order(-abs(get(paste0(g, "_logFC"))))]
  head(sub$gene_id, TOP_N)
})))
cat(sprintf("   union of top-%d DEGs across genotypes: %d genes\n", TOP_N, length(top_ids)))

mt <- m[gene_id %in% top_ids]
lfc_mat <- as.matrix(mt[, ..lfc_cols]); rownames(lfc_mat) <- mt$label_safe <- ifelse(
  !is.na(mt$acronym) & mt$acronym != "", mt$acronym, mt$gene_name)
colnames(lfc_mat) <- GENOTYPES
lfc_clip <- pmax(pmin(lfc_mat, LFC_CLIP), -LFC_CLIP)

# Significance stars matrix aligned to lfc_mat.
star_mat <- matrix("", nrow(lfc_mat), ncol(lfc_mat), dimnames = dimnames(lfc_mat))
for (j in seq_along(GENOTYPES)) {
  fdr <- mt[[fdr_cols[j]]]; sig <- mt[[sig_cols[j]]]
  star_mat[, j] <- ifelse(sig == "NS", "",
                   ifelse(fdr < 1e-3, "***", ifelse(fdr < 1e-2, "**", "*")))
}

show_rows_dd <- nrow(lfc_clip) <= 80
ph_dd <- pheatmap(lfc_clip,
  color = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100),
  breaks = seq(-LFC_CLIP, LFC_CLIP, length.out = 101),
  cluster_rows = TRUE, cluster_cols = TRUE, clustering_method = "ward.D2",
  display_numbers = if (show_rows_dd) star_mat else FALSE,
  number_color = "black", fontsize_number = 6,
  show_rownames = show_rows_dd, fontsize_row = 6, fontsize_col = 10,
  angle_col = 45, border_color = NA, na_col = "grey95",
  main = sprintf("Top DEGs (union of top %d / genotype) - log2FC vs WT", TOP_N),
  silent = TRUE)
hh_dd <- max(7, min(26, nrow(lfc_clip) * 0.12 + 3))
save_pheatmap(ph_dd, "Summary_Heatmap_topDEGs", width = 8, height = hh_dd)
fwrite(cbind(mt[, .(gene_id, gene_name, acronym, geneProduct)],
             as.data.table(lfc_mat)),
       file.path(SUM_DIR, "Summary_Heatmap_topDEGs.tsv"), sep = "\t")

# ===========================================================================
# 3. Curated symbiosis-panel heatmap (ported from J002, category-ordered)
# ===========================================================================
cat("== 3. Curated symbiosis-panel heatmap ==\n")
# category -> c(display_label = locus_tag). Locus tags verified present in the
# J003 annotation. Genes absent from the *filtered* set show as grey (NA).
key_genes <- list(
  "Leghemoglobin" = c(
    Lb2="MtrunA17Chr5g0427351", Lb1="MtrunA17Chr5g0435611",
    Lb5="MtrunA17Chr5g0417631", Lb7="MtrunA17Chr5g0435981",
    Lb9="MtrunA17Chr5g0435991", Lb11="MtrunA17Chr7g0270491"),
  "N-fixation / Symbiosis" = c(
    NAD1="MtrunA17Chr7g0222521", LYM1="MtrunA17Chr3g0114621",
    SymCRK="MtrunA17Chr3g0119041", NFS2="MtrunA17Chr8g0361201",
    NIN="MtrunA17Chr5g0448621", NIP="MtrunA17Chr1g0147631",
    CHIT5a="MtrunA17Chr1g0149821", CHIT5b="MtrunA17Chr4g0065511"),
  "NCR peptides" = c(
    NCR051="MtrunA17Chr5g0433121", NCR031="MtrunA17Chr5g0436091",
    NCR539="MtrunA17Chr7g0225561", NCR252="MtrunA17Chr2g0302921",
    NCR117="MtrunA17Chr5g0436571", NCR660="MtrunA17Chr5g0429331",
    NCR533="MtrunA17Chr7g0225531", NCR038="MtrunA17Chr5g0436051",
    NCR739="MtrunA17Chr8g0357471", NCR265="MtrunA17Chr5g0418951",
    NCR658="MtrunA17Chr8g0351581"),
  "Defense / Immunity" = c(
    `PR10-1`="MtrunA17Chr2g0295001", `PR10-4`="MtrunA17Chr2g0295141",
    PR3="MtrunA17Chr4g0026291", PR4="MtrunA17Chr5g0398561",
    CAPE12="MtrunA17Chr5g0411421", LRK="MtrunA17Chr2g0285861",
    DOX="MtrunA17Chr6g0451261"),
  "Senescence" = c(
    CP8="MtrunA17Chr3g0119321", `CP5/3`="MtrunA17Chr3g0119311",
    CP6="MtrunA17Chr4g0040881", CP2="MtrunA17Chr3g0119301",
    CHITINASE="MtrunA17Chr4g0027141", PAP="MtrunA17Chr3g0117691",
    VPE="MtrunA17Chr1g0151801"),
  "Redox" = c(
    `TRX S1`="MtrunA17Chr2g0317721", GRX9="MtrunA17Chr1g0195971")
)

labels   <- unlist(lapply(key_genes, names), use.names = FALSE)
locus    <- unlist(key_genes, use.names = FALSE)
category <- rep(names(key_genes), lengths(key_genes))

# Look up each locus tag by gene_name; build logFC + star matrices in fixed order.
idx <- match(locus, m$gene_name)
found <- !is.na(idx)
cat(sprintf("   panel genes present in filtered data: %d / %d\n", sum(found), length(locus)))
if (any(!found))
  cat("   (absent, shown grey):", paste(labels[!found], collapse = ", "), "\n")

cur_lfc <- matrix(NA_real_, length(labels), length(GENOTYPES),
                  dimnames = list(labels, GENOTYPES))
cur_star <- matrix("", length(labels), length(GENOTYPES), dimnames = dimnames(cur_lfc))
for (i in which(found)) {
  row <- m[idx[i]]
  for (j in seq_along(GENOTYPES)) {
    cur_lfc[i, j] <- row[[lfc_cols[j]]]
    if (row[[sig_cols[j]]] != "NS") {
      f <- row[[fdr_cols[j]]]
      cur_star[i, j] <- if (f < 1e-3) "***" else if (f < 1e-2) "**" else "*"
    }
  }
}
cur_clip <- pmax(pmin(cur_lfc, LFC_CLIP), -LFC_CLIP)

row_anno  <- data.frame(Category = factor(category, levels = names(key_genes)),
                        row.names = labels)
cat_pal   <- setNames(colorRampPalette(brewer.pal(8, "Set2"))(length(key_genes)),
                      names(key_genes))
gaps_rows <- head(cumsum(lengths(key_genes)), -1)

ph_cur <- pheatmap(cur_clip,
  color = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100),
  breaks = seq(-LFC_CLIP, LFC_CLIP, length.out = 101),
  cluster_rows = FALSE, cluster_cols = FALSE, gaps_row = gaps_rows,
  annotation_row = row_anno, annotation_colors = list(Category = cat_pal),
  annotation_names_row = FALSE,
  display_numbers = cur_star, number_color = "black", fontsize_number = 7,
  cellwidth = 26, cellheight = 13, fontsize_row = 8, fontsize_col = 10,
  angle_col = 45, border_color = "grey80", na_col = "grey92",
  legend_breaks = c(-4, -2, 0, 2, 4), legend_labels = c("-4","-2","0","+2","+4"),
  main = "Curated symbiosis panel - log2FC vs WT", silent = TRUE)
save_pheatmap(ph_cur, "Summary_Heatmap_curated", width = 7, height = 11)
fwrite(data.table(label = labels, category = category, locus_tag = locus,
                  as.data.table(cur_lfc)),
       file.path(SUM_DIR, "Summary_Heatmap_curated.tsv"), sep = "\t")

cat("\n== DONE ==\n")
cat("Outputs under:", SUM_DIR, "\n")
