#!/usr/bin/env Rscript
# =============================================================================
# 06.GO_KEGG_enrichment.R
# GO + KEGG over-representation analysis (ORA) of the 954 genes DOWN-regulated
# in BOTH A86485 and dnf2 (the "shared symbiotic-failure repression signature"
# identified by the A86485-vs-dnf2 Venn).
# -----------------------------------------------------------------------------
# DESIGN RATIONALE
#   * BACKGROUND = the 25,950 genes that were TESTED for DE (passed filterByExpr
#     in 00.edgeR_unified.R), NOT the whole genome. A gene that was never tested
#     could never have been called a DEG, so including it in the background would
#     inflate enrichment significance. This is the statistically correct universe.
#   * GO: direct (Blast2GO) annotations -> clusterProfiler::enricher, run
#     separately for BP / MF / CC. Term names + ontology come from GO.db.
#   * KEGG: KoFamScan assigns KEGG-Orthology (KO) numbers to proteins. We take the
#     single BEST significant KO per gene (the '*'-flagged, score>=threshold hit
#     with the highest score), map KO -> pathway via the live KEGG REST API, and
#     run pathway-level ORA at the GENE level (so the universe and the member
#     lists stay gene-based and interpretable). Global "overview" maps
#     (ko01100 Metabolic pathways, etc.) are excluded -- they are too generic.
#   * ID JOIN. Annotation files key genes as  MtrunA17_Chr3g0105081  (underscore);
#     our DE gene_name is  MtrunA17Chr3g0105081  (no underscore). The transform
#     is `MtrunA17` -> `MtrunA17_`, verified 1:1 and lossless. Hard assertions
#     guard it (this is the exact mismatch that caused an earlier 0/41 bug).
#
# COVERAGE CAVEAT (printed at runtime, important for interpretation)
#   Of the 954 DOWN genes, ~622 carry GO and only ~125 carry a KO. Many hallmark
#   nodule genes (NCR peptides, nodule-specific cysteine-rich genes, some
#   transporters) are functionally uncharacterized and therefore GO/KEGG-invisible.
#   ORA here surfaces the metabolic / defense / housekeeping arms of the response;
#   absence of a term is NOT absence of biology for that unannotated fraction.
#
# OUTPUT (99.Summary/Enrichment_A86485_dnf2_sharedDOWN/)
#   GO_enrichment_all.tsv        full GO ORA table (BP+MF+CC), member genes named
#   KEGG_enrichment.tsv          full KEGG pathway ORA table, member genes named
#   coverage_summary.tsv         annotation coverage of set & universe
#   gene2GO_used.tsv, gene2KO_used.tsv   the mappings actually used (provenance)
#   GO_dotplot_{BP,MF,CC}.*, GO_barplot_combined.*, GO_cnet_BP.*
#   KEGG_dotplot.*, KEGG_barplot.*, KEGG_cnet.*    (PDF + PNG + TIFF @600 dpi)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(clusterProfiler)
  library(GO.db)
  library(AnnotationDbi)
  library(enrichplot)
  library(ggplot2)
})

# ----------------------------- paths & config -------------------------------
BASE  <- "../05.DE_analysis"
SUMM  <- file.path(BASE, "99.Summary")
GLIST <- file.path(SUMM, "Venn_A86485_vs_dnf2", "GeneList_A86485_vs_dnf2.tsv")
MAST  <- file.path(SUMM, "MasterTable_with_DE.tsv")
ANNO  <- "${ANNOTATION_DIR}"
TOPGO <- file.path(ANNO, "sources/blast2go/MtrunA17r5.0-ANR-EGN-r1.9.b2g.topgo")
KOFAM <- file.path(ANNO, "annotation/MtrunA17r5.0-ANR-EGN-r1.9.prot.KoFamScan.txt")
OUT   <- file.path(SUMM, "Enrichment_A86485_dnf2_sharedDOWN")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

IMG_FORMATS <- c("pdf", "png", "tiff")
IMG_DPI     <- 900
FDR_CUT     <- 0.05          # significance threshold for reporting / plotting
GO_MIN <- 5;  GO_MAX <- 500  # gene-set size bounds for GO ORA
KG_MIN <- 3;  KG_MAX <- 500  # gene-set size bounds for KEGG ORA
TOP_N  <- 20                 # top terms to show in plots

# KEGG global/overview maps to drop (too generic for ORA)
KEGG_GLOBAL <- c("map01100","map01110","map01120","map01200","map01210",
                 "map01212","map01230","map01220","map01232","map01250","map01240")

# ID transforms between DE gene_name space and annotation-file space
to_annot <- function(x) sub("^MtrunA17", "MtrunA17_", x)   # MtrunA17Chr.. -> MtrunA17_Chr..

# Save a ggplot (enrichplot dotplot/cnet are ggplots too) to all formats.
save_plot <- function(p, base, width = 8, height = 7) {
  for (fmt in IMG_FORMATS) {
    fp <- file.path(OUT, paste0(base, ".", fmt))
    a  <- list(filename = fp, plot = p, width = width, height = height,
               dpi = IMG_DPI, device = fmt, limitsize = FALSE)
    if (fmt == "tiff") a$compression <- "lzw"
    do.call(ggsave, a)
  }
}

theme_pub <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

# =============================================================================
# 1. Load the query set (954 shared-DOWN) and the tested-gene universe
# =============================================================================
cat("== 1. Loading gene sets ==\n")
gl   <- fread(GLIST)
down <- gl[direction == "DOWN" & region == "A86485_dnf2_shared"]
stopifnot(nrow(down) == 954L)                     # known truth from the Venn

mast <- fread(MAST)
universe_name <- unique(mast$gene_name)
stopifnot(length(universe_name) == 25950L)        # filterByExpr-kept genes

sig_ids  <- to_annot(down$gene_name)              # query in annotation-ID space
univ_ids <- to_annot(universe_name)               # background in annotation-ID space
stopifnot(!any(duplicated(sig_ids)), !any(duplicated(univ_ids)))
stopifnot(all(sig_ids %in% univ_ids))             # DEGs must be a subset of universe
cat(sprintf("   query = %d shared-DOWN genes | universe = %d tested genes\n",
            length(sig_ids), length(univ_ids)))

# annotation-id -> readable label (acronym if present, else gene_name) for outputs
ann_map <- unique(mast[, .(annot = to_annot(gene_name), gene_name, acronym, geneProduct)])
ann_map[, label := fifelse(!is.na(acronym) & acronym != "", acronym, gene_name)]
lab_of <- setNames(ann_map$label, ann_map$annot)
# translate a "/"-separated clusterProfiler geneID string into readable labels
readable_ids <- function(s) vapply(strsplit(s, "/"), function(v)
  paste(ifelse(is.na(lab_of[v]), v, lab_of[v]), collapse = "/"), character(1))

# =============================================================================
# 2. GO ORA  (BP / MF / CC) via clusterProfiler::enricher
# =============================================================================
cat("== 2. GO enrichment ==\n")
tg <- fread(TOPGO, header = FALSE, sep = "\t", col.names = c("gene", "go"))
g2go <- tg[, .(GO = unlist(strsplit(go, ";"))), by = gene][GO != ""]
cat(sprintf("   b2g.topgo: %d genes, %d gene-GO pairs\n",
            uniqueN(g2go$gene), nrow(g2go)))

# attach ontology + human-readable term name from GO.db (drop obsolete GO ids)
all_go <- unique(g2go$GO)
go_info <- suppressWarnings(AnnotationDbi::select(
  GO.db, keys = all_go, columns = c("ONTOLOGY", "TERM"), keytype = "GOID"))
setDT(go_info)
go_info <- go_info[!is.na(ONTOLOGY)]
n_obsolete <- length(all_go) - uniqueN(go_info$GOID)
cat(sprintf("   GO ids resolved in GO.db: %d (%d obsolete/absent dropped)\n",
            uniqueN(go_info$GOID), n_obsolete))
g2go <- merge(g2go, go_info, by.x = "GO", by.y = "GOID")

run_go <- function(ont) {
  sub <- g2go[ONTOLOGY == ont]
  t2g <- sub[, .(term = GO, gene = gene)]
  t2n <- unique(sub[, .(term = GO, name = TERM)])
  enricher(gene = sig_ids, universe = univ_ids,
           TERM2GENE = t2g, TERM2NAME = t2n,
           pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "BH",
           minGSSize = GO_MIN, maxGSSize = GO_MAX)
}
go_res <- lapply(c("BP", "MF", "CC"), run_go)
names(go_res) <- c("BP", "MF", "CC")

# combined table with ontology + readable member genes
go_tab <- rbindlist(lapply(names(go_res), function(o) {
  r <- go_res[[o]]; if (is.null(r) || nrow(as.data.frame(r)) == 0) return(NULL)
  dt <- as.data.table(as.data.frame(r)); dt[, ONTOLOGY := o]; dt
}), fill = TRUE)
if (nrow(go_tab)) {
  go_tab[, geneSymbols := readable_ids(geneID)]
  setorder(go_tab, p.adjust)
  setcolorder(go_tab, c("ONTOLOGY","ID","Description","GeneRatio","BgRatio",
                        "pvalue","p.adjust","qvalue","Count","geneSymbols","geneID"))
  fwrite(go_tab, file.path(OUT, "GO_enrichment_all.tsv"), sep = "\t")
}
for (o in names(go_res)) {
  n_sig <- if (is.null(go_res[[o]])) 0 else sum(as.data.frame(go_res[[o]])$p.adjust < FDR_CUT)
  cat(sprintf("   %s: %d terms tested, %d significant (FDR<%.2f)\n",
              o, if (is.null(go_res[[o]])) 0 else nrow(as.data.frame(go_res[[o]])),
              n_sig, FDR_CUT))
}

# =============================================================================
# 3. KEGG ORA via KoFamScan (best significant KO per gene) + KEGG REST
# =============================================================================
cat("== 3. KEGG enrichment ==\n")
# significant KO assignments are the lines starting with '*' (score >= threshold)
ktab <- fread(cmd = sprintf("grep -E '^\\*' %s", shQuote(KOFAM)),
              header = FALSE, sep = "\n", quote = "")
parts <- tstrsplit(trimws(sub("^\\*", "", ktab$V1)), "\\s+")     # gene KO thr score eval def...
kdt   <- data.table(gene = parts[[1]], KO = parts[[2]],
                    score = as.numeric(parts[[4]]))
# BEST significant KO per gene = highest score
setorder(kdt, gene, -score)
g2ko <- kdt[, .SD[1], by = gene][, .(gene, KO)]
cat(sprintf("   KoFamScan: %d genes with a best significant KO\n", nrow(g2ko)))

# KO -> pathway and pathway names, fetched live from KEGG REST
fetch_rest <- function(url) {
  x <- fread(cmd = sprintf("curl -s --max-time 30 '%s'", url),
             header = FALSE, sep = "\t", quote = "")
  if (nrow(x) == 0) stop("KEGG REST returned no data: ", url)
  x
}
link <- fetch_rest("https://rest.kegg.jp/link/pathway/ko")
setnames(link, c("ko", "path"))
link <- link[grepl("^path:map", path)]                           # keep map##### (drop dup ko#####)
link[, KO := sub("^ko:", "", ko)][, pathway := sub("^path:", "", path)]
pname <- fetch_rest("https://rest.kegg.jp/list/pathway")
setnames(pname, c("pathway", "pname"))

# gene -> pathway (compose), drop global maps
g2path <- merge(g2ko, link[, .(KO, pathway)], by = "KO", allow.cartesian = TRUE)
g2path <- g2path[!pathway %in% KEGG_GLOBAL]
t2g_k <- unique(g2path[, .(term = pathway, gene = gene)])
t2n_k <- unique(pname[, .(term = pathway, name = pname)])
cat(sprintf("   KO->pathway: %d genes map to %d non-global pathways\n",
            uniqueN(g2path$gene), uniqueN(g2path$pathway)))

kegg_res <- enricher(gene = sig_ids, universe = univ_ids,
                     TERM2GENE = t2g_k, TERM2NAME = t2n_k,
                     pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "BH",
                     minGSSize = KG_MIN, maxGSSize = KG_MAX)
kegg_tab <- if (is.null(kegg_res)) data.table() else as.data.table(as.data.frame(kegg_res))
if (nrow(kegg_tab)) {
  kegg_tab[, geneSymbols := readable_ids(geneID)]
  setorder(kegg_tab, p.adjust)
  fwrite(kegg_tab, file.path(OUT, "KEGG_enrichment.tsv"), sep = "\t")
  cat(sprintf("   KEGG: %d pathways tested, %d significant (FDR<%.2f)\n",
              nrow(kegg_tab), sum(kegg_tab$p.adjust < FDR_CUT), FDR_CUT))
}

# =============================================================================
# 4. Coverage summary + provenance of the mappings used
# =============================================================================
cov <- data.table(
  set        = c("query_shared_DOWN", "tested_universe"),
  n_genes    = c(length(sig_ids), length(univ_ids)),
  n_with_GO  = c(sum(sig_ids %in% g2go$gene),  sum(univ_ids %in% g2go$gene)),
  n_with_KO  = c(sum(sig_ids %in% g2ko$gene),  sum(univ_ids %in% g2ko$gene)))
cov[, pct_GO := round(100 * n_with_GO / n_genes, 1)]
cov[, pct_KO := round(100 * n_with_KO / n_genes, 1)]
fwrite(cov, file.path(OUT, "coverage_summary.tsv"), sep = "\t")
print(cov)

# the gene->GO and gene->KO mappings restricted to the query set (provenance)
fwrite(merge(data.table(gene = sig_ids), g2go, by = "gene")[
         , .(gene, GO, ONTOLOGY, TERM)], file.path(OUT, "gene2GO_used.tsv"), sep = "\t")
fwrite(merge(data.table(gene = sig_ids), g2ko, by = "gene"),
       file.path(OUT, "gene2KO_used.tsv"), sep = "\t")

# =============================================================================
# 5. Figures
# =============================================================================
cat("== 5. Figures ==\n")
# enrichplot dotplots (one per GO ontology + KEGG), only when there are hits
dot_safe <- function(res, base, title) {
  if (is.null(res)) return(invisible())
  df <- as.data.frame(res); if (sum(df$p.adjust < FDR_CUT) == 0) {
    cat(sprintf("   (no FDR<%.2f hits; skipping %s)\n", FDR_CUT, base)); return(invisible()) }
  p <- dotplot(res, showCategory = TOP_N, x = "GeneRatio") +
    ggtitle(title) + theme_pub
  save_plot(p, base, width = 8.5, height = 8)
}
dot_safe(go_res$BP, "GO_dotplot_BP", "GO Biological Process - shared-DOWN")
dot_safe(go_res$MF, "GO_dotplot_MF", "GO Molecular Function - shared-DOWN")
dot_safe(go_res$CC, "GO_dotplot_CC", "GO Cellular Component - shared-DOWN")
dot_safe(kegg_res, "KEGG_dotplot",  "KEGG pathways - shared-DOWN")

# combined GO barplot: top terms per ontology, -log10 FDR, coloured by ontology
if (nrow(go_tab)) {
  bar <- go_tab[p.adjust < FDR_CUT]
  if (nrow(bar)) {
    bar <- bar[, head(.SD[order(p.adjust)], 10), by = ONTOLOGY]
    bar[, Description := factor(Description, levels = rev(unique(Description[order(ONTOLOGY, p.adjust)])))]
    g <- ggplot(bar, aes(-log10(p.adjust), Description, fill = ONTOLOGY)) +
      geom_col() + facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free_y") +
      scale_fill_manual(values = c(BP = "#1B9E77", MF = "#D95F02", CC = "#7570B3")) +
      labs(x = expression(-log[10]~FDR), y = NULL,
           title = "GO enrichment (top 10/ontology) - A86485 ∩ dnf2 shared-DOWN") +
      geom_vline(xintercept = -log10(FDR_CUT), linetype = 2, colour = "grey50") +
      theme_pub + theme(legend.position = "none")
    save_plot(g, "GO_barplot_combined", width = 9, height = 10)
  }
}

# cnetplot (gene-term network) with readable gene labels; defensive
cnet_safe <- function(res, base, title, n = 8) {
  if (is.null(res)) return(invisible())
  df <- as.data.frame(res); if (sum(df$p.adjust < FDR_CUT) == 0) return(invisible())
  r2 <- res
  r2@result$geneID <- readable_ids(r2@result$geneID)   # relabel members to acronym/name
  tryCatch({
    p <- cnetplot(r2, showCategory = n, node_label = "all",
                  cex_label_gene = 0.5) + ggtitle(title)
    save_plot(p, base, width = 11, height = 9)
  }, error = function(e) cat(sprintf("   (cnet %s skipped: %s)\n", base, conditionMessage(e))))
}
cnet_safe(go_res$BP, "GO_cnet_BP", "GO BP gene-concept network - shared-DOWN")
cnet_safe(kegg_res, "KEGG_cnet",  "KEGG gene-concept network - shared-DOWN")

cat("== DONE ==\n")
