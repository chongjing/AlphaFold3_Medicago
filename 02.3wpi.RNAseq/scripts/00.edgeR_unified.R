#!/usr/bin/env Rscript
# ============================================================================
# 00.edgeR_unified.R
# Differential expression: each genotype (3 reps) vs WT (3 reps), M. truncatula
# ----------------------------------------------------------------------------
# DESIGN RATIONALE (why this differs from the older per-pairwise template)
#
#   * UNIFIED FIT. One DGEList over ALL 21 samples with a ~0+group design;
#     each genotype-vs-WT contrast is then extracted from the single fit.
#     - Dispersion is estimated from WITHIN-group residuals across all groups,
#       giving 21-7 = 14 residual d.f. (vs only 4 d.f. in a 6-sample pairwise
#       fit) -> far more stable dispersion / QL shrinkage.
#     - Because the design is ~0+group, between-group mean differences do NOT
#       leak into the dispersion estimate; including a very divergent genotype
#       (e.g. dnf2) does not weaken the milder contrasts, provided each group's
#       own replicates are consistent.
#     - One common gene universe (single filterByExpr on the full design) makes
#       cross-comparison summaries (Venn/heatmap/merged table) apples-to-apples.
#
#   * glmTreat (FC-aware testing). Significance is tested AGAINST the fold-change
#     threshold (|log2FC| > LFC_THRESH) rather than against 0 followed by a post
#     hoc logFC filter. glmTreat controls the error rate for the fold-change
#     claim itself and yields fewer, higher-confidence DEGs (edgeR User's Guide).
#
#   * Annotation enrichment. Curated acronym / geneProduct are merged in for
#     readable plot labels; genes without a curated entry keep their locus tag.
#
#   * Outputs are written as TSV (never comma-quoted CSV) because geneProduct
#     descriptions contain commas that would corrupt a bare CSV.
#
# INPUTS  (produced by 05.build_expression_tables.R)
#   <MATRIX_DIR>/counts_matrix.tsv     raw htseq counts, genes x 21 samples
#   <MATRIX_DIR>/TPM_matrix.tsv        StringTie TPM   (carried into result tables)
#   <MATRIX_DIR>/FPKM_matrix.tsv       StringTie FPKM  (carried into result tables)
#   <MATRIX_DIR>/gene_annotation.tsv   gene_id, gene_name, chr, strand, start, end, length
#   <FUNC_ANNOT>                       locus_tag, acronym, geneProduct (curated)
#
# OUTPUTS
#   00.QC/                             global QC over all 21 samples
#       MDS, PCA(logCPM), sample-correlation heatmap, library-size barplot
#   <genotype>_vs_WT/                  one folder per contrast
#       01.BCV / 02.QL_dispersion / 03.Volcano / 04.MA / 05.Heatmap /
#       06.Pvalue_hist  + DE_results_full.tsv / DE_results_significant.tsv /
#       DE_summary.tsv
#   99.Summary/
#       DEG_counts_barplot.*           UP/DOWN per genotype
#       All_comparisons_DE_summary.tsv stacked summary
#       MasterTable_with_DE.tsv        annotation + per-genotype logFC/FDR/sig
#
# USAGE
#   Rscript 00.edgeR_unified.R                 # all default contrasts
#   Rscript 00.edgeR_unified.R A256 dnf2       # only the named genotypes
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
})

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
BASE_DIR   <- ".."
MATRIX_DIR <- file.path(BASE_DIR, "03.Mapping", "matrix")
OUT_DIR    <- file.path(BASE_DIR, "05.DE_analysis")
FUNC_ANNOT <- "${ANNOTATION_DIR}/annotation/MtrunA17r5.0-ANR-EGN-r1.9.gene_annotation.tsv"

REFERENCE   <- "WT"                 # baseline group
LFC_THRESH  <- 1.0                  # |log2FC| threshold tested by glmTreat
FDR_CUTOFF  <- 0.05                 # significance cutoff on the glmTreat FDR
IMG_FORMATS <- c("pdf", "png", "tiff")
IMG_DPI     <- 600

# Replicate layout: every sample column is "<group>_<rep>".
N_REP <- 3

# ---------------------------------------------------------------------------
# Plot helpers
# ---------------------------------------------------------------------------
theme_pub <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      axis.title    = element_text(size = rel(1.0), face = "bold"),
      axis.text     = element_text(size = rel(0.85), color = "black"),
      plot.title    = element_text(size = rel(1.1), face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = rel(0.9), hjust = 0.5, color = "grey40"),
      axis.line     = element_line(linewidth = 0.5, color = "black"),
      axis.ticks    = element_line(linewidth = 0.5, color = "black"),
      legend.title  = element_text(size = rel(0.9), face = "bold"),
      legend.text   = element_text(size = rel(0.85)),
      legend.key.size = unit(0.8, "lines"),
      panel.border  = element_blank(),
      panel.background = element_blank(),
      plot.margin   = margin(6, 6, 6, 6)
    )
}

# Save a ggplot object to all requested raster/vector formats.
# (compression is only a valid device arg for tiff, so pass it conditionally.)
save_ggplot <- function(p, dir, base, width = 7, height = 6) {
  for (fmt in IMG_FORMATS) {
    fp <- file.path(dir, paste0(base, ".", fmt))
    a <- list(filename = fp, plot = p, width = width, height = height,
              dpi = IMG_DPI, device = fmt)
    if (fmt == "tiff") a$compression <- "lzw"
    do.call(ggsave, a)
  }
}

# Open a device of the given format (helper shared by the savers below).
open_dev <- function(fp, fmt, width, height) {
  if (fmt == "pdf")  pdf(fp, width = width, height = height)
  if (fmt == "png")  png(fp, width = width, height = height, units = "in", res = IMG_DPI)
  if (fmt == "tiff") tiff(fp, width = width, height = height, units = "in",
                          res = IMG_DPI, compression = "lzw")
}

# Save a base-graphics expression (edgeR plotBCV/plotQLDisp) to all formats.
save_base <- function(expr, dir, base, width = 6, height = 5) {
  expr <- substitute(expr)
  for (fmt in IMG_FORMATS) {
    open_dev(file.path(dir, paste0(base, ".", fmt)), fmt, width, height)
    eval(expr, envir = parent.frame())
    dev.off()
  }
}

# Save a pheatmap object to all formats. pheatmap(silent=TRUE) returns the
# heatmap but does NOT draw it, so we must grid.draw its $gtable explicitly;
# eval'ing the call alone yields blank pages.
save_pheatmap <- function(ph, dir, base, width = 8, height = 7) {
  for (fmt in IMG_FORMATS) {
    open_dev(file.path(dir, paste0(base, ".", fmt)), fmt, width, height)
    grid::grid.newpage()
    grid::grid.draw(ph$gtable)
    dev.off()
  }
}

# ===========================================================================
# 1. Load matrices and annotation; assert gene_id integrity throughout
# ===========================================================================
cat("== 1. Loading matrices ==\n")
counts <- fread(file.path(MATRIX_DIR, "counts_matrix.tsv"))
tpm    <- fread(file.path(MATRIX_DIR, "TPM_matrix.tsv"))
fpkm   <- fread(file.path(MATRIX_DIR, "FPKM_matrix.tsv"))
annot  <- fread(file.path(MATRIX_DIR, "gene_annotation.tsv"))

# All four tables must describe the SAME genes in the SAME order. We do not
# trust row order: assert identical gene_id vectors and abort otherwise.
stopifnot(identical(counts$gene_id, tpm$gene_id),
          identical(counts$gene_id, fpkm$gene_id),
          identical(counts$gene_id, annot$gene_id))
cat(sprintf("   %d genes x %d samples\n", nrow(counts), ncol(counts) - 1L))

# Merge curated functional annotation (acronym, geneProduct). The locus_tag in
# the curated file carries an extra underscore ("MtrunA17_Chr..") that our
# gene_name lacks; normalize before joining. Left join => all genes retained.
fa <- fread(FUNC_ANNOT, sep = "\t", quote = "")
fa[, gene_name := gsub("MtrunA17_", "MtrunA17", locus_tag)]
fa <- unique(fa, by = "gene_name")
annot <- merge(annot, fa[, .(gene_name, acronym, geneProduct)],
               by = "gene_name", all.x = TRUE, sort = FALSE)
# Re-key to counts order (merge may reorder) and re-assert.
annot <- annot[match(counts$gene_id, gene_id)]
stopifnot(identical(annot$gene_id, counts$gene_id))
# A readable label: curated acronym when present, else the locus tag.
annot[, label := ifelse(!is.na(acronym) & acronym != "", acronym, gene_name)]
cat(sprintf("   curated acronyms matched: %d\n",
            sum(!is.na(annot$acronym) & annot$acronym != "")))

# ---------------------------------------------------------------------------
# Build the count matrix (genes x samples) with gene_id rownames.
# ---------------------------------------------------------------------------
sample_cols <- setdiff(colnames(counts), "gene_id")
mat <- as.matrix(counts[, ..sample_cols])
rownames(mat) <- counts$gene_id

# Group = sample name minus the trailing "_<rep>".
group_of <- sub("_[0-9]+$", "", sample_cols)
groups   <- unique(group_of)
stopifnot(REFERENCE %in% groups)
treatments <- setdiff(groups, REFERENCE)

# Optional CLI restriction to a subset of genotypes.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  bad <- setdiff(args, treatments)
  if (length(bad)) stop("Unknown genotype(s): ", paste(bad, collapse = ", "),
                        "\n  available: ", paste(treatments, collapse = ", "))
  treatments <- args
}
cat(sprintf("   reference: %s | contrasts: %s\n",
            REFERENCE, paste(treatments, collapse = ", ")))

# ===========================================================================
# 2. Unified DGEList, filtering, TMM normalization
# ===========================================================================
cat("== 2. DGEList + filterByExpr + TMM (all samples) ==\n")
group <- factor(group_of, levels = c(REFERENCE, sort(treatments)))
y <- DGEList(counts = mat, group = group)

design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

# Single filter on the full design => one common gene universe for every contrast.
keep <- filterByExpr(y, design)
cat(sprintf("   genes: %d -> %d kept (%.1f%%)\n",
            nrow(y), sum(keep), 100 * mean(keep)))
y <- y[keep, , keep.lib.sizes = FALSE]
y <- normLibSizes(y)

# ===========================================================================
# 3. Dispersion + QL fit (once, shared by all contrasts)
# ===========================================================================
cat("== 3. estimateDisp + glmQLFit (robust) ==\n")
y   <- estimateDisp(y, design, robust = TRUE)
fit <- glmQLFit(y, design, robust = TRUE)
cat(sprintf("   common dispersion = %.4f | BCV = %.4f | residual df = %d\n",
            y$common.dispersion, sqrt(y$common.dispersion),
            ncol(mat) - ncol(design)))

# logCPM (normalized) once; reused for QC and per-gene expression columns.
logcpm <- cpm(y, log = TRUE)

# ===========================================================================
# 4. Global QC over all 21 samples (run once)
# ===========================================================================
cat("== 4. Global QC ==\n")
qc_dir <- file.path(OUT_DIR, "00.QC")
dir.create(qc_dir, showWarnings = FALSE, recursive = TRUE)

grp_levels <- levels(group)
grp_pal <- setNames(colorRampPalette(brewer.pal(8, "Dark2"))(length(grp_levels)),
                    grp_levels)

# --- 4a. MDS (edgeR leading-logFC) ---
mds <- plotMDS(y, plot = FALSE)
mds_df <- data.frame(Dim1 = mds$x, Dim2 = mds$y,
                     Sample = colnames(y), Group = as.character(group))
p_mds <- ggplot(mds_df, aes(Dim1, Dim2, color = Group, label = Sample)) +
  geom_point(size = 3.5, alpha = 0.9) +
  geom_text_repel(size = 2.6, show.legend = FALSE, max.overlaps = Inf) +
  scale_color_manual(values = grp_pal) +
  labs(title = "MDS (leading logFC) - all samples",
       x = sprintf("Dim 1 (%.1f%%)", mds$var.explained[1] * 100),
       y = sprintf("Dim 2 (%.1f%%)", mds$var.explained[2] * 100)) +
  theme_pub()
save_ggplot(p_mds, qc_dir, "01.MDS_all_samples", width = 7, height = 6)

# --- 4b. PCA on logCPM ---
pca <- prcomp(t(logcpm), scale. = FALSE)
ve  <- 100 * pca$sdev^2 / sum(pca$sdev^2)
pca_df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                     Sample = colnames(y), Group = as.character(group))
p_pca <- ggplot(pca_df, aes(PC1, PC2, color = Group, label = Sample)) +
  geom_point(size = 3.5, alpha = 0.9) +
  geom_text_repel(size = 2.6, show.legend = FALSE, max.overlaps = Inf) +
  scale_color_manual(values = grp_pal) +
  labs(title = "PCA (logCPM) - all samples",
       x = sprintf("PC1 (%.1f%%)", ve[1]), y = sprintf("PC2 (%.1f%%)", ve[2])) +
  theme_pub()
save_ggplot(p_pca, qc_dir, "02.PCA_logCPM", width = 7, height = 6)

# --- 4c. Sample-sample correlation heatmap (Spearman on logCPM) ---
cor_mat <- cor(logcpm, method = "spearman")
ann_col <- data.frame(Group = as.character(group)); rownames(ann_col) <- colnames(y)
ph_cor <- pheatmap(cor_mat, annotation_col = ann_col, annotation_row = ann_col,
                   annotation_colors = list(Group = grp_pal),
                   color = colorRampPalette(brewer.pal(9, "Blues"))(100),
                   display_numbers = FALSE,
                   main = "Sample correlation (Spearman, logCPM)", silent = TRUE)
save_pheatmap(ph_cor, qc_dir, "03.Sample_correlation", width = 8, height = 7)

# --- 4d. Library size barplot ---
libdf <- data.frame(Sample = colnames(y),
                    LibSize = y$samples$lib.size / 1e6,
                    Group = as.character(group))
libdf$Sample <- factor(libdf$Sample, levels = libdf$Sample)
p_lib <- ggplot(libdf, aes(Sample, LibSize, fill = Group)) +
  geom_col() + scale_fill_manual(values = grp_pal) +
  labs(title = "Library size (post-filter)", y = "Million reads", x = NULL) +
  theme_pub() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_ggplot(p_lib, qc_dir, "04.Library_sizes", width = 8, height = 5)

# ===========================================================================
# 5. Per-contrast DE (glmTreat against LFC_THRESH)
# ===========================================================================
# Static annotation block carried into every result table, aligned to the
# filtered gene set by gene_id (NOT by position).
annot_keep <- annot[match(rownames(y), gene_id)]
stopifnot(identical(annot_keep$gene_id, rownames(y)))

summary_rows <- list()
master_de    <- data.table(gene_id = rownames(y))

for (g in treatments) {
  cat(sprintf("== 5. Contrast: %s vs %s ==\n", g, REFERENCE))
  cdir <- file.path(OUT_DIR, paste0(g, "_vs_", REFERENCE))
  dir.create(cdir, showWarnings = FALSE, recursive = TRUE)

  # contrast = genotype - reference (vector over design columns)
  ctr <- setNames(rep(0, ncol(design)), colnames(design))
  ctr[g] <- 1; ctr[REFERENCE] <- -1

  # glmTreat: test |log2FC| > LFC_THRESH directly.
  res <- glmTreat(fit, contrast = ctr, lfc = LFC_THRESH)
  tt  <- topTags(res, n = Inf, sort.by = "none")$table
  tt  <- tt[rownames(y), ]              # enforce gene order == y

  # Per-sample logCPM for this contrast's samples (ref + this genotype).
  this_samples <- colnames(y)[group %in% c(REFERENCE, g)]
  cpm_block <- logcpm[, this_samples, drop = FALSE]
  colnames(cpm_block) <- paste0(this_samples, "_logCPM")

  de <- data.table(
    gene_id    = rownames(y),
    gene_name  = annot_keep$gene_name,
    label      = annot_keep$label,
    acronym    = annot_keep$acronym,
    geneProduct= annot_keep$geneProduct,
    chr        = annot_keep$chr,
    strand     = annot_keep$strand,
    start      = annot_keep$start,
    end        = annot_keep$end,
    logFC      = tt$logFC,
    logCPM     = tt$logCPM,
    PValue     = tt$PValue,
    FDR        = tt$FDR
  )
  de <- cbind(de, as.data.table(cpm_block))

  # Significance: glmTreat FDR already accounts for the FC threshold, so the
  # call is simply FDR < cutoff; direction from the sign of logFC.
  de[, change := "NS"]
  de[FDR < FDR_CUTOFF & logFC > 0,  change := "UP"]
  de[FDR < FDR_CUTOFF & logFC < 0,  change := "DOWN"]
  n_up <- sum(de$change == "UP"); n_dn <- sum(de$change == "DOWN")
  cat(sprintf("   UP=%d  DOWN=%d  (FDR<%.2f, glmTreat lfc=%.1f)\n",
              n_up, n_dn, FDR_CUTOFF, LFC_THRESH))

  # ---- write tables (TSV) ----
  fwrite(de, file.path(cdir, "DE_results_full.tsv"), sep = "\t")
  fwrite(de[change != "NS"], file.path(cdir, "DE_results_significant.tsv"), sep = "\t")

  # ---- diagnostics: BCV + QL dispersion (shared model, drawn per folder) ----
  save_base(plotBCV(y, main = sprintf("BCV: %s vs %s", g, REFERENCE)),
            cdir, "01.BCV", width = 6, height = 5)
  save_base(plotQLDisp(fit, main = sprintf("QL dispersion: %s vs %s", g, REFERENCE)),
            cdir, "02.QL_dispersion", width = 6, height = 5)

  # ---- volcano ----
  cols <- c(UP = "#B2182B", DOWN = "#2166AC", NS = "grey70")
  top <- de[change != "NS"][order(FDR)][seq_len(min(20, .N))]
  p_v <- ggplot(de, aes(logFC, -log10(FDR), color = change)) +
    geom_point(size = 1.2, alpha = 0.7, shape = 16) +
    scale_color_manual(values = cols,
                       breaks = c("UP", "DOWN", "NS"),
                       labels = c(sprintf("UP (%d)", n_up),
                                  sprintf("DOWN (%d)", n_dn), "NS")) +
    geom_vline(xintercept = c(-LFC_THRESH, LFC_THRESH), linetype = "dashed",
               color = "grey40", linewidth = 0.4) +
    geom_hline(yintercept = -log10(FDR_CUTOFF), linetype = "dashed",
               color = "grey40", linewidth = 0.4) +
    geom_text_repel(data = top, aes(label = label), size = 2.5,
                    max.overlaps = 20, show.legend = FALSE,
                    box.padding = 0.5, segment.color = "grey50") +
    labs(title = sprintf("Volcano: %s vs %s", g, REFERENCE),
         subtitle = sprintf("UP %d | DOWN %d | FDR<%.2f, glmTreat |log2FC|>%.1f",
                            n_up, n_dn, FDR_CUTOFF, LFC_THRESH),
         x = expression(log[2]~"fold change"), y = expression(-log[10]~"FDR"),
         color = NULL) +
    theme_pub() +
    theme(legend.position = "inside", legend.position.inside = c(0.99, 0.99),
          legend.justification = c(1, 1))
  save_ggplot(p_v, cdir, "03.Volcano", width = 7, height = 6)

  # ---- MA ----
  p_ma <- ggplot(de, aes(logCPM, logFC, color = change)) +
    geom_point(size = 1.0, alpha = 0.6, shape = 16) +
    scale_color_manual(values = cols) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.3) +
    geom_hline(yintercept = c(-LFC_THRESH, LFC_THRESH), linetype = "dashed",
               color = "grey40", linewidth = 0.4) +
    labs(title = sprintf("MA: %s vs %s", g, REFERENCE),
         x = expression(log[2]~"CPM"), y = expression(log[2]~"fold change"),
         color = NULL) +
    theme_pub() +
    theme(legend.position = "inside", legend.position.inside = c(0.99, 0.99),
          legend.justification = c(1, 1))
  save_ggplot(p_ma, cdir, "04.MA", width = 7, height = 6)

  # ---- heatmap of significant DEGs (row-scaled logCPM) ----
  sig_idx <- which(de$change != "NS")
  if (length(sig_idx) >= 2) {
    hm <- cpm_block[sig_idx, , drop = FALSE]
    colnames(hm) <- sub("_logCPM$", "", colnames(hm))
    hm_scaled <- t(scale(t(hm)))
    hm_scaled <- hm_scaled[stats::complete.cases(hm_scaled), , drop = FALSE]
    show_rows <- nrow(hm_scaled) <= 60
    hann <- data.frame(Group = ifelse(grepl(paste0("^", REFERENCE), colnames(hm_scaled)),
                                      REFERENCE, g))
    rownames(hann) <- colnames(hm_scaled)
    hcol <- grp_pal[c(REFERENCE, g)]
    hh <- max(6, min(22, nrow(hm_scaled) / 30 + 4))
    ph_deg <- pheatmap(hm_scaled,
               color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
               annotation_col = hann, annotation_colors = list(Group = hcol),
               show_rownames = show_rows,
               labels_row = if (show_rows) annot_keep$label[sig_idx] else NULL,
               clustering_method = "ward.D2",
               main = sprintf("%s vs %s - %d DEGs", g, REFERENCE, nrow(hm_scaled)),
               silent = TRUE)
    save_pheatmap(ph_deg, cdir, "05.Heatmap_DEGs", width = 8, height = hh)
  } else {
    cat("   (<2 DEGs; skipping heatmap)\n")
  }

  # ---- p-value histogram (diagnostic) ----
  p_h <- ggplot(de, aes(PValue)) +
    geom_histogram(bins = 50, fill = "grey60", color = "white", linewidth = 0.3) +
    labs(title = sprintf("P-value distribution: %s vs %s", g, REFERENCE),
         subtitle = "Flat under null; peak near 0 = true signal",
         x = "P-value", y = "Count") +
    theme_pub()
  save_ggplot(p_h, cdir, "06.Pvalue_hist", width = 6, height = 5)

  # ---- per-contrast summary + accumulate cross-comparison structures ----
  fwrite(data.table(Comparison = sprintf("%s_vs_%s", g, REFERENCE),
                    Genes_tested = nrow(de), UP = n_up, DOWN = n_dn,
                    Total_DEG = n_up + n_dn,
                    LFC_threshold = LFC_THRESH, FDR_cutoff = FDR_CUTOFF),
         file.path(cdir, "DE_summary.tsv"), sep = "\t")

  summary_rows[[g]] <- data.table(Comparison = sprintf("%s_vs_%s", g, REFERENCE),
                                  Genotype = g, UP = n_up, DOWN = n_dn,
                                  Total_DEG = n_up + n_dn)
  master_de[, (paste0(g, "_logFC")) := de$logFC]
  master_de[, (paste0(g, "_FDR"))   := de$FDR]
  master_de[, (paste0(g, "_sig"))   := de$change]
}

# ===========================================================================
# 6. Cross-comparison summary
# ===========================================================================
cat("== 6. Cross-comparison summary ==\n")
sum_dir <- file.path(OUT_DIR, "99.Summary")
dir.create(sum_dir, showWarnings = FALSE, recursive = TRUE)

all_sum <- rbindlist(summary_rows)
fwrite(all_sum, file.path(sum_dir, "All_comparisons_DE_summary.tsv"), sep = "\t")

# Merged master DE table: static annotation + per-genotype logFC/FDR/sig.
master_out <- cbind(annot_keep[, .(gene_id, gene_name, acronym, geneProduct,
                                   chr, strand, start, end)],
                    master_de[, -1])
fwrite(master_out, file.path(sum_dir, "MasterTable_with_DE.tsv"), sep = "\t")

# DEG-count barplot (UP positive, DOWN negative) ordered by total.
bar_df <- rbind(
  data.table(Genotype = all_sum$Genotype, n = all_sum$UP,   dir = "UP"),
  data.table(Genotype = all_sum$Genotype, n = -all_sum$DOWN, dir = "DOWN")
)
ord <- all_sum[order(-Total_DEG)]$Genotype
bar_df$Genotype <- factor(bar_df$Genotype, levels = ord)
p_bar <- ggplot(bar_df, aes(Genotype, n, fill = dir)) +
  geom_col(width = 0.7) +
  geom_text(data = all_sum, aes(Genotype, UP, label = UP),
            inherit.aes = FALSE, vjust = -0.3, size = 3) +
  geom_text(data = all_sum, aes(Genotype, -DOWN, label = DOWN),
            inherit.aes = FALSE, vjust = 1.2, size = 3) +
  scale_fill_manual(values = c(UP = "#B2182B", DOWN = "#2166AC")) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
  labs(title = sprintf("DEGs vs %s (FDR<%.2f, glmTreat |log2FC|>%.1f)",
                       REFERENCE, FDR_CUTOFF, LFC_THRESH),
       x = NULL, y = "DOWN          # DEGs          UP", fill = NULL) +
  theme_pub() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_ggplot(p_bar, sum_dir, "DEG_counts_barplot", width = 8, height = 6)

cat("\n== DONE ==\n")
print(all_sum)
cat(sprintf("Outputs under: %s\n", OUT_DIR))
