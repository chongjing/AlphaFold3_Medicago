#!/usr/bin/env python3
# ============================================================================
# 03.venn_diagrams.py
# Area-proportional 3-way Venn diagrams of DEG overlaps, built from the
# MasterTable_with_DE.tsv produced by 00.edgeR_unified.R.
#
# The genotype names are combinatorial:
#     A86256 = A86 + A256   (double)
#     A86485 = A86 + A485   (double)
# so the biologically meaningful 3-way comparisons are each double mutant
# against its two constituent singles. Each tests whether the double mutant's
# transcriptome response is the union of the singles (additive), a subset
# (redundancy/epistasis), or contains novel genes (synergy).
#
# For each triple we draw UP and DOWN separately (a gene up in one set and
# down in another should not be merged), and write the 7 region gene lists +
# a summary CSV.
#
# matplotlib_venn lives in the base anaconda python (not the R envs), so run
# with that interpreter:
#   python3 03.venn_diagrams.py
# ============================================================================

import csv
import os
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib_venn import venn3, venn3_circles

BASE = ".."
MASTER = os.path.join(BASE, "05.DE_analysis", "99.Summary", "MasterTable_with_DE.tsv")
OUT_DIR = os.path.join(BASE, "05.DE_analysis", "99.Summary", "Venn")
os.makedirs(OUT_DIR, exist_ok=True)

DPI = 600
FORMATS = ["pdf", "png", "tiff"]

# Combinatorial triples: (set_A_single, set_B_single, double_mutant)
TRIPLES = [
    ("A86", "A256", "A86256"),
    ("A86", "A485", "A86485"),
]

# Colorblind-safe fills for the three circles.
CIRCLE_COLORS = ("#4DBBD5", "#E64B35", "#00A087")


def load_sets(path):
    """Return {genotype: {'UP': set(gene_id), 'DOWN': set(gene_id)}} from the
    *_sig columns of the master table."""
    with open(path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    genotypes = [c[:-4] for c in rows[0].keys() if c.endswith("_sig")]
    sets = {g: {"UP": set(), "DOWN": set()} for g in genotypes}
    for r in rows:
        for g in genotypes:
            s = r.get(f"{g}_sig", "NS")
            if s in ("UP", "DOWN"):
                sets[g][s].add(r["gene_id"])
    return sets, len(rows)


def venn_regions(A, B, C):
    """7 disjoint regions keyed by membership code (matplotlib_venn order)."""
    return {
        "100": A - B - C,
        "010": B - A - C,
        "001": C - A - B,
        "110": (A & B) - C,
        "101": (A & C) - B,
        "011": (B & C) - A,
        "111": A & B & C,
    }


def draw(sets, triple, direction, out_dir):
    a, b, c = triple
    A = sets[a][direction]; B = sets[b][direction]; C = sets[c][direction]
    regions = venn_regions(A, B, C)
    subset_sizes = {k: len(v) for k, v in regions.items()}

    fig, ax = plt.subplots(figsize=(6, 6))
    v = venn3(subsets=tuple(subset_sizes[k] for k in
              ["100", "010", "110", "001", "101", "011", "111"]),
              set_labels=(a, b, c), ax=ax,
              set_colors=CIRCLE_COLORS, alpha=0.55)
    venn3_circles(subsets=tuple(subset_sizes[k] for k in
                  ["100", "010", "110", "001", "101", "011", "111"]),
                  ax=ax, linewidth=1.0, color="grey")
    # Bold the subset count labels for readability.
    for sid in ["100", "010", "001", "110", "101", "011", "111"]:
        lbl = v.get_label_by_id(sid)
        if lbl is not None:
            lbl.set_fontsize(11)
    for lbl in v.set_labels:
        if lbl is not None:
            lbl.set_fontsize(13)
            lbl.set_fontweight("bold")

    color = "#B2182B" if direction == "UP" else "#2166AC"
    ax.set_title(f"{direction} DEGs vs WT\n{a} ∩ {b} ∩ {c}",
                 fontsize=13, fontweight="bold", color=color)

    base = os.path.join(out_dir, f"Venn_{a}_{b}_{c}_{direction}")
    for fmt in FORMATS:
        fig.savefig(f"{base}.{fmt}", dpi=DPI, bbox_inches="tight",
                    **({"pil_kwargs": {"compression": "tiff_lzw"}} if fmt == "tiff" else {}))
    plt.close(fig)

    # Write each region's gene list.
    for code, genes in regions.items():
        tag = {"100": a + "_only", "010": b + "_only", "001": c + "_only",
               "110": f"{a}_{b}_shared", "101": f"{a}_{c}_shared",
               "011": f"{b}_{c}_shared", "111": "all_three"}[code]
        with open(f"{base}__{tag}.txt", "w") as fh:
            fh.write("\n".join(sorted(genes)) + ("\n" if genes else ""))
    return subset_sizes


def main():
    sets, n_genes = load_sets(MASTER)
    print(f"Loaded master table: {n_genes} genes")
    summary = []
    for triple in TRIPLES:
        a, b, c = triple
        for direction in ("UP", "DOWN"):
            sizes = draw(sets, triple, direction, OUT_DIR)
            print(f"  {a}/{b}/{c} {direction}: "
                  + ", ".join(f"{k}={v}" for k, v in sizes.items()))
            row = {"triple": f"{a}/{b}/{c}", "direction": direction}
            row.update({
                f"{a}_only": sizes["100"], f"{b}_only": sizes["010"],
                f"{c}_only": sizes["001"], f"{a}&{b}": sizes["110"],
                f"{a}&{c}": sizes["101"], f"{b}&{c}": sizes["011"],
                "all_three": sizes["111"],
            })
            summary.append(row)

    # One tidy CSV (union of keys across both triples).
    keys = []
    for r in summary:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(os.path.join(OUT_DIR, "Venn_summary.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        for r in summary:
            w.writerow(r)
    print(f"\nDone. Outputs under: {OUT_DIR}")


if __name__ == "__main__":
    main()
