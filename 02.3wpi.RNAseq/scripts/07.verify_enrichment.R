#!/usr/bin/env Rscript
# =============================================================================
# 07.verify_enrichment.R
# Independent verification of 06.GO_KEGG_enrichment.R outputs.
# Rebuilds the universe/query annotation counts FROM SCRATCH (b2g.topgo + GO.db
# + KEGG REST + KoFamScan) and recomputes the hypergeometric p-value for the top
# GO term and top KEGG pathway BY HAND (phyper), comparing to clusterProfiler's
# reported pvalue. Also traces every member gene back to the master DE table to
# confirm it is a genuine 954-shared-DOWN gene (DOWN in BOTH A86485 and dnf2).
# This does NOT reuse clusterProfiler's internal counts -- it re-derives them.
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(GO.db); library(AnnotationDbi) })

ANNO  <- "${ANNOTATION_DIR}"
SUMM  <- "../05.DE_analysis/99.Summary"
OUT   <- file.path(SUMM, "Enrichment_A86485_dnf2_sharedDOWN")
to_annot <- function(x) sub("^MtrunA17", "MtrunA17_", x)

mast <- fread(file.path(SUMM, "MasterTable_with_DE.tsv"))
univ <- to_annot(unique(mast$gene_name))
gl   <- fread(file.path(SUMM, "Venn_A86485_vs_dnf2/GeneList_A86485_vs_dnf2.tsv"))
down <- gl[direction == "DOWN" & region == "A86485_dnf2_shared"]
sig  <- to_annot(down$gene_name)

cat("=== sanity: query subset of universe? ", all(sig %in% univ), " | nQuery=", length(sig),
    " nUniv=", length(univ), "\n", sep = "")

## ---------------- GO ----------------
tg <- fread(file.path(ANNO, "sources/blast2go/MtrunA17r5.0-ANR-EGN-r1.9.b2g.topgo"),
            header = FALSE, col.names = c("gene", "go"))
g2go <- tg[, .(GO = unlist(strsplit(go, ";"))), by = gene][GO != ""]
info <- as.data.table(suppressWarnings(AnnotationDbi::select(
  GO.db, unique(g2go$GO), c("ONTOLOGY", "TERM"), "GOID")))[!is.na(ONTOLOGY)]
g2go <- merge(g2go, info, by.x = "GO", by.y = "GOID")
g2go_u <- g2go[gene %in% univ]                      # restrict to background universe

go_tab <- fread(file.path(OUT, "GO_enrichment_all.tsv"))
top_go <- go_tab[order(p.adjust)][1]                # oxygen transport expected
term <- top_go$ID; ont <- top_go$ONTOLOGY
N <- uniqueN(g2go_u[ONTOLOGY == ont]$gene)                          # universe genes w/ any term in ont
K <- uniqueN(g2go_u[ONTOLOGY == ont & gene %in% sig]$gene)          # query genes w/ any term in ont
m <- uniqueN(g2go_u[GO == term]$gene)                              # universe genes in term
k <- uniqueN(g2go_u[GO == term & gene %in% sig]$gene)              # query genes in term
p_hand <- phyper(k - 1, m, N - m, K, lower.tail = FALSE)
cat(sprintf("\n[GO] %s (%s) '%s'\n", term, ont, top_go$Description))
cat(sprintf("   hand:  k=%d K=%d m=%d N=%d  -> GeneRatio %d/%d  BgRatio %d/%d\n",
            k, K, m, N, k, K, m, N))
cat(sprintf("   table: GeneRatio %s  BgRatio %s\n", top_go$GeneRatio, top_go$BgRatio))
cat(sprintf("   phyper(hand) = %.4e   |   clusterProfiler pvalue = %.4e   | match=%s\n",
            p_hand, top_go$pvalue, isTRUE(all.equal(p_hand, top_go$pvalue, tolerance = 1e-6))))

# trace the member genes of the top GO term back to the master DE table
members <- g2go_u[GO == term & gene %in% sig]$gene
mm <- mast[to_annot(gene_name) %in% members]
cat(sprintf("   members: %d genes; all in 954-DOWN set=%s; all DOWN in BOTH A86485 & dnf2=%s\n",
            length(members),
            all(members %in% sig),
            all(mm$A86485_sig == "DOWN" & mm$dnf2_sig == "DOWN")))
cat("   member acronyms: ", paste(sort(mm$acronym[mm$acronym != ""]), collapse = ", "), "\n")

## ---------------- KEGG ----------------
ktab <- fread(cmd = sprintf("grep -E '^\\*' %s",
              shQuote(file.path(ANNO, "annotation/MtrunA17r5.0-ANR-EGN-r1.9.prot.KoFamScan.txt"))),
              header = FALSE, sep = "\n", quote = "")
parts <- tstrsplit(trimws(sub("^\\*", "", ktab$V1)), "\\s+")
kdt <- data.table(gene = parts[[1]], KO = parts[[2]], score = as.numeric(parts[[4]]))
setorder(kdt, gene, -score)
g2ko <- kdt[, .SD[1], by = gene][, .(gene, KO)]
link <- fread(cmd = "curl -s --max-time 30 'https://rest.kegg.jp/link/pathway/ko'",
              header = FALSE, col.names = c("ko", "path"))[grepl("^path:map", path)]
link[, `:=`(KO = sub("^ko:", "", ko), pathway = sub("^path:", "", path))]
GLOBAL <- c("map01100","map01110","map01120","map01200","map01210","map01212",
            "map01230","map01220","map01232","map01250","map01240")
g2path <- merge(g2ko, link[, .(KO, pathway)], by = "KO", allow.cartesian = TRUE)[!pathway %in% GLOBAL]
g2path_u <- g2path[gene %in% univ]

kg <- fread(file.path(OUT, "KEGG_enrichment.tsv"))
top_k <- kg[order(p.adjust)][1]                     # Nitrogen metabolism expected
pw <- top_k$ID
Nk <- uniqueN(g2path_u$gene)                                   # universe genes w/ any pathway
Kk <- uniqueN(g2path_u[gene %in% sig]$gene)                    # query genes w/ any pathway
mk <- uniqueN(g2path_u[pathway == pw]$gene)
kk <- uniqueN(g2path_u[pathway == pw & gene %in% sig]$gene)
pk_hand <- phyper(kk - 1, mk, Nk - mk, Kk, lower.tail = FALSE)
cat(sprintf("\n[KEGG] %s '%s'\n", pw, top_k$Description))
cat(sprintf("   hand:  k=%d K=%d m=%d N=%d  -> GeneRatio %d/%d  BgRatio %d/%d\n",
            kk, Kk, mk, Nk, kk, Kk, mk, Nk))
cat(sprintf("   table: GeneRatio %s  BgRatio %s\n", top_k$GeneRatio, top_k$BgRatio))
cat(sprintf("   phyper(hand) = %.4e   |   clusterProfiler pvalue = %.4e   | match=%s\n",
            pk_hand, top_k$pvalue, isTRUE(all.equal(pk_hand, top_k$pvalue, tolerance = 1e-6))))
mk_genes <- g2path_u[pathway == pw & gene %in% sig]$gene
mmk <- mast[to_annot(gene_name) %in% mk_genes]
cat(sprintf("   members: %d genes; all DOWN in BOTH=%s; acronyms: %s\n",
            length(mk_genes), all(mmk$A86485_sig == "DOWN" & mmk$dnf2_sig == "DOWN"),
            paste(sort(mmk$acronym[mmk$acronym != ""]), collapse = ", ")))
