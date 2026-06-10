#!/usr/bin/env python3
# ============================================================================
# 06.venn_pairwise.py
# Three pairwise 2-way Venn diagrams of DEG overlaps, built from the
# MasterTable_with_DE.tsv produced by 00.edgeR_unified.R.
#
# Comparisons (all vs WT):
#   1. A86  vs A485   — two single mutants
#   2. A86  vs dnf2   — single mutant vs Fix- control
#   3. A485 vs dnf2   — single mutant vs Fix- control
#
# Follows the exact pattern of 04.venn_A86485_vs_dnf2.py:
#   - area-proportional 2-way Venn (matplotlib_venn)
#   - UP and DOWN drawn separately
#   - annotated gene list (every gene in any region)
#   - summary CSV per comparison
#
# Run with the base anaconda python (matplotlib_venn lives there, not the R envs):
#   python3 06.venn_pairwise.py
# ============================================================================

import csv
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn2_circles

BASE = ".."
MASTER = os.path.join(BASE, "05.DE_analysis", "99.Summary", "MasterTable_with_DE.tsv")
DPI = 600
FORMATS = ["pdf", "png", "tiff"]
COLORS = ("#E64B35", "#3C5488")   # A/warm, B/cool  (same as 04.venn)

# The three pairwise comparisons.
PAIRS = [
    ("A86",  "A485"),
    ("A86",  "dnf2"),
    ("A485", "dnf2"),
]


def load_master(path):
    """Load MasterTable_with_DE.tsv as list of dicts."""
    with open(path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return rows


def sig_set(rows, genotype, direction):
    """Return set of gene_ids called as `direction` (UP/DOWN) for `genotype`."""
    return {r["gene_id"] for r in rows if r.get(f"{genotype}_sig") == direction}


def draw_venn(A_set, B_set, A_name, B_name, direction, out_dir):
    """Draw a 2-way Venn, save in all formats, return (only_A, only_B, shared)."""
    only_A = len(A_set - B_set)
    only_B = len(B_set - A_set)
    shared = len(A_set & B_set)

    fig, ax = plt.subplots(figsize=(6.5, 6))
    v = venn2(subsets=(only_A, only_B, shared), set_labels=(A_name, B_name), ax=ax,
              set_colors=COLORS, alpha=0.6)
    venn2_circles(subsets=(only_A, only_B, shared), ax=ax,
                  linewidth=1.2, color="grey")

    # Plain count inside each region.
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
    ax.set_title(f"{direction} DEGs vs WT — {A_name} ∩ {B_name}",
                 fontsize=14, fontweight="bold", color=color, pad=16)

    # Percentage-of-set caption beneath the diagram.
    if A_set and B_set:
        ax.text(0.5, -0.06,
                f"shared = {shared}  "
                f"({100*shared/len(A_set):.0f}% of {A_name}, "
                f"{100*shared/len(B_set):.0f}% of {B_name})",
                transform=ax.transAxes, ha="center", va="top",
                fontsize=11, color="#333333")

    base = os.path.join(out_dir, f"Venn_{A_name}_vs_{B_name}_{direction}")
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
        return "shared"
    if in_a:
        return "A_only"
    if in_b:
        return "B_only"
    return None


def process_pair(rows, A_name, B_name):
    """Run one pairwise comparison: Venn diagrams + annotated gene list + summary."""
    out_dir = os.path.join(BASE, "05.DE_analysis", "99.Summary",
                           f"Venn_{A_name}_vs_{B_name}")
    os.makedirs(out_dir, exist_ok=True)

    print(f"\n{'='*60}")
    print(f"  {A_name} vs {B_name}")
    print(f"{'='*60}")

    summary = []
    for direction in ("UP", "DOWN"):
        A_set = sig_set(rows, A_name, direction)
        B_set = sig_set(rows, B_name, direction)
        oA, oB, sh = draw_venn(A_set, B_set, A_name, B_name, direction, out_dir)
        print(f"  {direction}: {A_name}_only={oA}  {B_name}_only={oB}  shared={sh}")
        summary.append({
            "direction": direction,
            f"{A_name}_only": oA, f"{B_name}_only": oB, "shared": sh,
            f"{A_name}_total": oA + sh, f"{B_name}_total": oB + sh,
        })

    # ---- annotated gene list: every gene appearing in ANY region ----
    ann_cols = ["gene_id", "gene_name", "acronym", "geneProduct",
                "chr", "strand", "start", "end"]
    de_cols = [f"{A_name}_logFC", f"{A_name}_FDR", f"{A_name}_sig",
               f"{B_name}_logFC", f"{B_name}_FDR", f"{B_name}_sig"]

    out_rows = []
    for r in rows:
        for direction in ("UP", "DOWN"):
            reg = region_of(r.get(f"{A_name}_sig"), r.get(f"{B_name}_sig"), direction)
            if reg is None:
                continue
            rec = {c: r.get(c, "") for c in ann_cols + de_cols}
            rec["direction"] = direction
            # Use A_name/B_name in region label for clarity.
            if reg == "shared":
                rec["region"] = f"{A_name}_{B_name}_shared"
            elif reg == "A_only":
                rec["region"] = f"{A_name}_only"
            else:
                rec["region"] = f"{B_name}_only"
            out_rows.append(rec)

    # Sort: direction, region, then strongest |B logFC| first for browsing.
    def keyfn(x):
        try:
            mag = -abs(float(x[f"{B_name}_logFC"]))
        except (ValueError, TypeError):
            mag = 0.0
        return (x["direction"], x["region"], mag)
    out_rows.sort(key=keyfn)

    list_path = os.path.join(out_dir, f"GeneList_{A_name}_vs_{B_name}.tsv")
    fieldnames = ["direction", "region"] + ann_cols + de_cols
    with open(list_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        w.writerows(out_rows)
    print(f"  annotated gene list ({len(out_rows)} rows) -> {os.path.basename(list_path)}")

    # ---- summary CSV ----
    with open(os.path.join(out_dir, f"Venn_{A_name}_vs_{B_name}_summary.csv"),
              "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(summary[0].keys()))
        w.writeheader()
        w.writerows(summary)

    return summary


def main():
    rows = load_master(MASTER)
    print(f"Loaded master table: {len(rows)} genes")

    all_summaries = {}
    for A_name, B_name in PAIRS:
        s = process_pair(rows, A_name, B_name)
        all_summaries[f"{A_name}_vs_{B_name}"] = s

    # ---- combined summary across all pairs ----
    combined_path = os.path.join(BASE, "05.DE_analysis", "99.Summary",
                                 "Venn_pairwise_combined_summary.csv")
    fieldnames = ["comparison", "direction",
                  "A_only", "B_only", "shared", "A_total", "B_total"]
    with open(combined_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        for pair_key, sums in all_summaries.items():
            A_name, B_name = pair_key.split("_vs_")
            for s in sums:
                row = {"comparison": pair_key, "direction": s["direction"],
                       "A_only": s[f"{A_name}_only"], "B_only": s[f"{B_name}_only"],
                       "shared": s["shared"],
                       "A_total": s[f"{A_name}_total"],
                       "B_total": s[f"{B_name}_total"]}
                w.writerow(row)
    print(f"\nCombined summary -> {os.path.basename(combined_path)}")
    print("\nDone.")


if __name__ == "__main__":
    main()
