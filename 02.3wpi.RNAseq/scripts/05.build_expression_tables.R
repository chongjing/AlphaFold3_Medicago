#!/usr/bin/env Rscript
# ============================================================================
# 05.build_expression_tables.R
#
# Purpose
# -------
# For each RNA-seq sample, merge the per-gene outputs produced during mapping
# into ONE "master" expression table, then combine all samples into wide
# gene x sample matrices for downstream DE analysis (DESeq2/edgeR, etc.).
#
# Inputs (located in each  <Mapping>/<sample>/  directory)
#   3.<sample>.count.tsv : htseq-count raw counts. Headerless, 2 columns:
#                          gene_id <tab> raw_count. The last 5 lines are
#                          htseq summary rows (__no_feature, __ambiguous, ...)
#                          and are dropped.
#   4.<sample>.FPKM.tsv  : StringTie -A gene abundance. Header present:
#                          "Gene ID, Gene Name, Reference, Strand, Start, End,
#                           Coverage, FPKM, TPM" (one row per gene).
#   t_data.ctab          : StringTie -B transcript table. Used only to obtain a
#                          true exonic transcript length per gene (col "length").
#
# Outputs
#   per sample : <sample>/<sample>.master_expression.tsv
#                (gene_id, gene_name, chr, strand, start, end, length,
#                 count, coverage, FPKM, TPM) -- one row per gene
#   combined   : <Mapping>/matrix/counts_matrix.tsv  (genes x samples, htseq counts)
#                <Mapping>/matrix/TPM_matrix.tsv     (genes x samples, StringTie TPM)
#                <Mapping>/matrix/FPKM_matrix.tsv    (genes x samples, StringTie FPKM)
#                <Mapping>/matrix/gene_annotation.tsv (gene_id + static annotation)
#
# IMPORTANT -- gene_id integrity
# ------------------------------
# Raw counts (htseq) and normalized values (StringTie) come from two DIFFERENT
# tools. They MUST be aligned strictly by gene_id, never by row order. This
# script therefore asserts, at every step, that the gene_id sets are identical
# and aborts (stop()) on any mismatch. A silent mis-alignment can never reach
# the output. Do NOT "fix" an abort by switching to positional binding.
#
# NOTE on the exon table (e_data.ctab): it is intentionally NOT used for
# gene-level counts. Its `rcount` counts every read OVERLAPPING each exon, so
# summing exons multi-counts multi-exon reads (verified: ~100x inflated vs
# htseq). htseq raw count is the correct gene-level count for DE analysis.
# ============================================================================

suppressPackageStartupMessages(library(data.table))

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
map_dir <- "../03.Mapping"
out_dir <- file.path(map_dir, "matrix")   # combined matrices go here

# Explicit sample list -> controls the column order of the matrices.
samples <- c(
  "A256_1",  "A256_2",  "A256_3",
  "A485_1",  "A485_2",  "A485_3",
  "A86_1",   "A86_2",   "A86_3",
  "A86256_1","A86256_2","A86256_3",
  "A86485_1","A86485_2","A86485_3",
  "dnf2_1",  "dnf2_2",  "dnf2_3",
  "WT_1",    "WT_2",    "WT_3"
)

# ---------------------------------------------------------------------------
# Helper: read one sample, merge its three files into a master table.
# Returns a data.table keyed by gene_id with all per-gene columns.
# ---------------------------------------------------------------------------
read_sample <- function(sample) {
  sdir <- file.path(map_dir, sample)

  ht_file <- file.path(sdir, paste0("3.", sample, ".count.tsv"))
  st_file <- file.path(sdir, paste0("4.", sample, ".FPKM.tsv"))
  td_file <- file.path(sdir, "t_data.ctab")
  for (f in c(ht_file, st_file, td_file)) {
    if (!file.exists(f)) stop("[", sample, "] missing input file: ", f)
  }

  # --- (1) htseq-count raw counts -----------------------------------------
  # Headerless: gene_id <tab> count. Drop the 5 trailing "__*" summary rows.
  ht <- fread(ht_file, header = FALSE, col.names = c("gene_id", "count"),
              colClasses = list(character = 1, integer = 2))
  ht <- ht[!grepl("^__", gene_id)]
  if (anyDuplicated(ht$gene_id))
    stop("[", sample, "] duplicate gene_id in htseq counts")

  # --- (2) StringTie -A gene abundance (FPKM / TPM / coords) ---------------
  st <- fread(st_file, header = TRUE)
  # Rename the StringTie headers (which contain spaces) to clean names.
  setnames(st,
           old = c("Gene ID","Gene Name","Reference","Strand",
                   "Start","End","Coverage","FPKM","TPM"),
           new = c("gene_id","gene_name","chr","strand",
                   "start","end","coverage","FPKM","TPM"))
  if (anyDuplicated(st$gene_id))
    stop("[", sample, "] duplicate gene_id in StringTie -A output")

  # --- (3) Gene length from t_data (longest transcript per gene) -----------
  # This annotation is 1 transcript per gene, but max() is used so the script
  # stays correct if a multi-transcript annotation is ever used.
  td <- fread(td_file, header = TRUE, select = c("gene_id", "length"))
  glen <- td[, .(length = max(length)), by = gene_id]

  # --- gene_id integrity check #1: htseq vs StringTie must be identical ----
  if (!setequal(ht$gene_id, st$gene_id)) {
    only_ht <- setdiff(ht$gene_id, st$gene_id)
    only_st <- setdiff(st$gene_id, ht$gene_id)
    stop("[", sample, "] gene_id mismatch htseq vs StringTie | ",
         "htseq-only=", length(only_ht), " StringTie-only=", length(only_st),
         " e.g. ", paste(head(c(only_ht, only_st), 3), collapse = ", "))
  }

  # --- merge, always BY gene_id (never by position) ------------------------
  # Inner joins; counts above already guarantee the sets are equal, so no rows
  # are dropped. The post-merge row-count assertion below makes that explicit.
  m <- merge(st, ht,   by = "gene_id", all = FALSE)
  m <- merge(m,  glen, by = "gene_id", all = FALSE)
  if (nrow(m) != nrow(ht))
    stop("[", sample, "] row count changed after merge (",
         nrow(ht), " -> ", nrow(m), ") -- gene_id join problem")

  # Final column order for the per-sample master table.
  setcolorder(m, c("gene_id","gene_name","chr","strand","start","end",
                   "length","count","coverage","FPKM","TPM"))
  setkey(m, gene_id)
  m
}

# ---------------------------------------------------------------------------
# Pass 1: build + write each per-sample master table.
# ---------------------------------------------------------------------------
cat("Building per-sample master tables...\n")
tables <- vector("list", length(samples))
names(tables) <- samples

for (s in samples) {
  m <- read_sample(s)
  out <- file.path(map_dir, s, paste0(s, ".master_expression.tsv"))
  fwrite(m, out, sep = "\t")
  cat(sprintf("  %-10s %6d genes -> %s\n", s, nrow(m), basename(out)))
  tables[[s]] <- m
}

# ---------------------------------------------------------------------------
# gene_id integrity check #2: ALL samples must share the same gene_id set.
# Use sample 1 as the reference order; every other sample must match exactly.
# ---------------------------------------------------------------------------
ref_ids <- tables[[1]]$gene_id
for (s in samples[-1]) {
  if (!setequal(tables[[s]]$gene_id, ref_ids))
    stop("[", s, "] gene_id set differs from reference sample '", samples[1],
         "' -- cannot build a consistent matrix")
}
cat(sprintf("\nAll %d samples share an identical set of %d gene_ids.\n",
            length(samples), length(ref_ids)))

# ---------------------------------------------------------------------------
# Pass 2: assemble wide matrices. Each column is realigned to `ref_ids` by
# gene_id (via keyed lookup), NOT by row order, then asserted NA-free.
# ---------------------------------------------------------------------------
extract_aligned <- function(dt, value_col) {
  v <- dt[.(ref_ids), on = "gene_id", get(value_col)]  # reorder to ref_ids
  if (anyNA(v))
    stop("alignment produced NA in column '", value_col,
         "' -- gene_id lookup failed")
  v
}

counts_mat <- data.table(gene_id = ref_ids)
tpm_mat    <- data.table(gene_id = ref_ids)
fpkm_mat   <- data.table(gene_id = ref_ids)
for (s in samples) {
  counts_mat[, (s) := extract_aligned(tables[[s]], "count")]
  tpm_mat[,    (s) := extract_aligned(tables[[s]], "TPM")]
  fpkm_mat[,   (s) := extract_aligned(tables[[s]], "FPKM")]
}

# Static per-gene annotation (identical across samples -> taken from reference).
annotation <- tables[[1]][, .(gene_id, gene_name, chr, strand, start, end, length)]

# ---------------------------------------------------------------------------
# Write combined outputs.
# ---------------------------------------------------------------------------
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
fwrite(counts_mat, file.path(out_dir, "counts_matrix.tsv"), sep = "\t")
fwrite(tpm_mat,    file.path(out_dir, "TPM_matrix.tsv"),    sep = "\t")
fwrite(fpkm_mat,   file.path(out_dir, "FPKM_matrix.tsv"),   sep = "\t")
fwrite(annotation, file.path(out_dir, "gene_annotation.tsv"), sep = "\t")

cat("\nCombined matrices written to:", out_dir, "\n")
cat(sprintf("  counts_matrix.tsv   %d genes x %d samples\n", nrow(counts_mat), length(samples)))
cat(sprintf("  TPM_matrix.tsv      %d genes x %d samples\n", nrow(tpm_mat),    length(samples)))
cat(sprintf("  FPKM_matrix.tsv     %d genes x %d samples\n", nrow(fpkm_mat),   length(samples)))
cat(sprintf("  gene_annotation.tsv %d genes\n", nrow(annotation)))
cat("\nDone.\n")
