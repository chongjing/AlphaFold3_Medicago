#!/usr/bin/env python3
# ============================================================================
# 04.venn_A86485_vs_dnf2.py
# Two-way DEG overlap: A86485 (double mutant) vs dnf2 (Fix- control mutant).
#
# Motivation: in the global PCA and the data-driven summary heatmap, A86485 and
# dnf2 cluster together, suggesting A86485 phenocopies the dnf2 senescence/Fix-
# response more than it resembles its A485 parent. This script quantifies that
# overlap directly (UP and DOWN separately) and exports an annotated gene list.
#
# Outputs (under 99.Summary/Venn_A86485_vs_dnf2/):
#   Venn_A86485_vs_dnf2_UP.{pdf,png,tiff}
#   Venn_A86485_vs_dnf2_DOWN.{pdf,png,tiff}
#   GeneList_A86485_vs_dnf2.tsv   annotated, every gene in any region, with
#                                 region/direction + logFC/FDR from BOTH genotypes
#   Venn_A86485_vs_dnf2_summary.csv
#
# Run with the base anaconda python (matplotlib_venn lives there, not the R envs):
#   python3 04.venn_A86485_vs_dnf2.py
# ============================================================================

import csv
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn2_circles

BASE = ".."
MASTER = os.path.join(BASE, "05.DE_analysis", "99.Summary", "MasterTable_with_DE.tsv")
OUT_DIR = os.path.join(BASE, "05.DE_analysis", "99.Summary", "Venn_A86485_vs_dnf2")
os.makedirs(OUT_DIR, exist_ok=True)

A, B = "A86485", "dnf2"           # the two genotypes to compare
DPI = 600
FORMATS = ["pdf", "png", "tiff"]
COLORS = ("#E64B35", "#3C5488")   # A86485 = warm, dnf2 = cool


def load_master(path):
    with open(path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return rows


def sig_set(rows, genotype, direction):
    return {r["gene_id"] for r in rows if r.get(f"{genotype}_sig") == direction}


def draw_venn(rows, direction):
    A_set = sig_set(rows, A, direction)
    B_set = sig_set(rows, B, direction)
    only_A = len(A_set - B_set)
    only_B = len(B_set - A_set)
    shared = len(A_set & B_set)

    fig, ax = plt.subplots(figsize=(6.5, 6))
    v = venn2(subsets=(only_A, only_B, shared), set_labels=(A, B), ax=ax,
              set_colors=COLORS, alpha=0.6)
    venn2_circles(subsets=(only_A, only_B, shared), ax=ax,
                  linewidth=1.2, color="grey")

    # Plain count inside each region; keep the intersection label to the bare
    # number so it never collides with the only-A crescent on asymmetric pairs.
    region_text = {"10": str(only_A), "01": str(only_B), "11": str(shared)}
    for sid, txt in region_text.items():
        lab = v.get_label_by_id(sid)
        if lab is not None:
            lab.set_text(txt)
            lab.set_fontsize(13)
            lab.set_fontweight("bold")
    for lab in v.set_labels:
        if lab is not None:
            lab.set_fontsize(14)
            lab.set_fontweight("bold")

    color = "#B2182B" if direction == "UP" else "#2166AC"
    ax.set_title(f"{direction} DEGs vs WT — {A} ∩ {B}",
                 fontsize=14, fontweight="bold", color=color, pad=16)
    # Percentage-of-set context as a caption beneath the diagram (avoids
    # cramming multi-line text into a small intersection lens).
    if A_set and B_set:
        ax.text(0.5, -0.06,
                f"shared = {shared}  "
                f"({100*shared/len(A_set):.0f}% of {A}, "
                f"{100*shared/len(B_set):.0f}% of {B})",
                transform=ax.transAxes, ha="center", va="top",
                fontsize=11, color="#333333")

    base = os.path.join(OUT_DIR, f"Venn_{A}_vs_{B}_{direction}")
    for fmt in FORMATS:
        kw = {"pil_kwargs": {"compression": "tiff_lzw"}} if fmt == "tiff" else {}
        fig.savefig(f"{base}.{fmt}", dpi=DPI, bbox_inches="tight", **kw)
    plt.close(fig)
    return only_A, only_B, shared


def region_of(a_sig, b_sig, direction):
    """Region label for a gene given its sig calls in A and B for `direction`."""
    in_a = a_sig == direction
    in_b = b_sig == direction
    if in_a and in_b:
        return f"{A}_{B}_shared"
    if in_a:
        return f"{A}_only"
    if in_b:
        return f"{B}_only"
    return None


def main():
    rows = load_master(MASTER)
    print(f"Loaded master table: {len(rows)} genes")

    summary = []
    for direction in ("UP", "DOWN"):
        oA, oB, sh = draw_venn(rows, direction)
        print(f"  {direction}: {A}_only={oA}  {B}_only={oB}  shared={sh}")
        summary.append({"direction": direction, f"{A}_only": oA,
                        f"{B}_only": oB, "shared": sh,
                        f"{A}_total": oA + sh, f"{B}_total": oB + sh})

    # ---- annotated gene list: every gene appearing in ANY region ----
    ann_cols = ["gene_id", "gene_name", "acronym", "geneProduct",
                "chr", "strand", "start", "end"]
    de_cols = [f"{A}_logFC", f"{A}_FDR", f"{A}_sig",
               f"{B}_logFC", f"{B}_FDR", f"{B}_sig"]
    out_rows = []
    for r in rows:
        for direction in ("UP", "DOWN"):
            reg = region_of(r.get(f"{A}_sig"), r.get(f"{B}_sig"), direction)
            if reg is None:
                continue
            rec = {c: r.get(c, "") for c in ann_cols + de_cols}
            rec["direction"] = direction
            rec["region"] = reg
            out_rows.append(rec)

    # sort: direction, region, then strongest |dnf2 logFC| first for browsing
    def keyfn(x):
        try:
            mag = -abs(float(x[f"{B}_logFC"]))
        except (ValueError, TypeError):
            mag = 0.0
        return (x["direction"], x["region"], mag)
    out_rows.sort(key=keyfn)

    list_path = os.path.join(OUT_DIR, f"GeneList_{A}_vs_{B}.tsv")
    fieldnames = ["direction", "region"] + ann_cols + de_cols
    with open(list_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        w.writerows(out_rows)
    print(f"  annotated gene list ({len(out_rows)} rows) -> {os.path.basename(list_path)}")

    with open(os.path.join(OUT_DIR, f"Venn_{A}_vs_{B}_summary.csv"),
              "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(summary[0].keys()))
        w.writeheader()
        w.writerows(summary)

    print(f"\nDone. Outputs under: {OUT_DIR}")


if __name__ == "__main__":
    main()
