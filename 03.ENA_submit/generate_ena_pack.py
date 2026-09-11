#!/usr/bin/env python3
"""Build the local ENA submission pack for Medicago DNF2-SRP nodule RNA-seq (J003).

Reads MD5.txt + os.stat (never hand-transcribes checksums/sizes), maps the 7
conditions to host-genotype x rhizobial-strain, and emits:
  - metadata.tsv                 (master table, 21 libraries)
  - sample_checklist_ERC000037.tsv
  - study_text.md
  - FASTQ symlinks in this dir (no copy)

No network. No submission. Study/SAMPLE accessions are empty until registration.
"""
from __future__ import annotations

import csv
import hashlib
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DATA = ROOT.parent / "X101SC26030664-Z01-J003"
RAW = DATA / "01.RawData"
MD5FILE = DATA / "MD5.txt"

STUDY_ALIAS = "Medicago_DNF2_SRP_RNAseq"
HOLD_DATE = "2026-11-01"
INSTRUMENT = "DNBSEQ-T7"
LIBRARY_SOURCE = "TRANSCRIPTOMIC"
LIBRARY_SELECTION = "Oligo-dT"
LIBRARY_STRATEGY = "RNA-Seq"
LAYOUT = "PAIRED"

# host genotype -> short label used in titles/attributes
HOST_LABEL = {
    "WT": "wild type",
    "dnf2": "dnf2 (DNF2 loss-of-function mutant)",
}

# rhizobial strain key -> human label
STRAIN_LABEL = {
    "WT": "wild type",
    "d-rsp86": "\u0394rsp86 deletion mutant",
    "d-rsp256": "\u0394rsp256 deletion mutant",
    "d-rsp485": "\u0394rsp485 deletion mutant",
    "d-rsp86 d-rsp256": "\u0394rsp86 \u0394rsp256 double deletion mutant",
    "d-rsp86 d-rsp485": "\u0394rsp86 \u0394rsp485 double deletion mutant",
}

# library filename prefix -> (host genotype key, rhizobial strain key)
CONDITIONS = {
    "WT": ("WT", "WT"),
    "dnf2": ("dnf2", "WT"),
    "A86": ("WT", "d-rsp86"),
    "A256": ("WT", "d-rsp256"),
    "A485": ("WT", "d-rsp485"),
    "A86256": ("WT", "d-rsp86 d-rsp256"),
    "A86485": ("WT", "d-rsp86 d-rsp485"),
}


def md5_of(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_md5file() -> dict[str, str]:
    """Return {filename: md5hex} for 01.RawData/*.fq.gz entries only."""
    out: dict[str, str] = {}
    for line in MD5FILE.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        parts = line.split()
        if len(parts) != 2:
            continue
        md5, rel = parts
        if "/01.RawData/" not in rel or not rel.endswith(".fq.gz"):
            continue
        out[rel.split("/")[-1]] = md5
    return out


def build_rows(md5map: dict[str, str]) -> list[dict]:
    rows = []
    # deterministic order: WT first, then dnf2, then the mutants, each rep 1..3
    order = ["WT", "dnf2", "A86", "A256", "A485", "A86256", "A86485"]
    for prefix in order:
        host_key, strain_key = CONDITIONS[prefix]
        for rep in (1, 2, 3):
            lib = f"{prefix}_{rep}"
            fq1 = f"{lib}_1.fq.gz"
            fq2 = f"{lib}_2.fq.gz"
            if fq1 not in md5map or fq2 not in md5map:
                raise SystemExit(f"missing MD5 for {fq1}/{fq2}")
            p1 = RAW / fq1
            p2 = RAW / fq2
            if not p1.exists() or not p2.exists():
                raise SystemExit(f"missing file {fq1} or {fq2}")
            rows.append({
                "library_name": lib,
                "sample_alias": f"MtrDNF2_{lib}",
                "host_genotype": host_key,
                "strain": strain_key,
                "replicate": str(rep),
                "instrument": INSTRUMENT,
                "library_source": LIBRARY_SOURCE,
                "library_selection": LIBRARY_SELECTION,
                "library_strategy": LIBRARY_STRATEGY,
                "layout": LAYOUT,
                "fq1": fq1, "fq1_md5": md5map[fq1], "fq1_bytes": str(p1.stat().st_size),
                "fq2": fq2, "fq2_md5": md5map[fq2], "fq2_bytes": str(p2.stat().st_size),
                "study_accession": "", "sample_accession": "",
            })
    return rows


FIELDS = [
    "library_name", "sample_alias", "host_genotype", "strain", "replicate",
    "instrument", "library_source", "library_selection", "library_strategy",
    "layout", "fq1", "fq1_md5", "fq1_bytes", "fq2", "fq2_md5", "fq2_bytes",
    "study_accession", "sample_accession",
]


def write_metadata(rows: list[dict]) -> None:
    with (ROOT / "metadata.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"Wrote metadata.tsv n={len(rows)}")


def write_checklist(rows: list[dict]) -> None:
    hdr = [
        "sample_alias", "tax_id", "scientific_name", "common_name",
        "sample_title", "sample_description", "collection date",
        "geographic location (country and/or sea)",
        "geographic location (region and locality)", "cultivar",
        "genotype", "organism part", "isolation_source", "infect",
        "collected_by",
    ]
    lines = ["#checklist_accession\tERC000037",
             "#" + "\t".join(hdr)]
    for r in rows:
        title, desc = sample_text(r)
        lines.append("\t".join([
            r["sample_alias"], "3880", "Medicago truncatula",
            "barrel medic", title, desc, "not collected",
            "United Kingdom", "not collected", "A17",
            HOST_LABEL[r["host_genotype"]], "root nodule", "root nodule",
            "Sinorhizobium meliloti" if r["strain"] != "WT" else "none",
            "Chongjing Xia",
        ]))
    (ROOT / "sample_checklist_ERC000037.tsv").write_text("\n".join(lines) + "\n")
    print(f"Wrote sample_checklist_ERC000037.tsv n={len(rows)}")


def sample_text(r: dict) -> tuple[str, str]:
    host = HOST_LABEL[r["host_genotype"]]
    strain = STRAIN_LABEL[r["strain"]]
    title = (f"Medicago truncatula A17 root nodule, {host}, inoculated with "
             f"Sinorhizobium meliloti 2011 {strain}, replicate {r['replicate']}")
    desc = (
        f"Medicago truncatula (ecotype A17) root nodule RNA-seq. Host genotype: "
        f"{host}. Plants germinated on nitrogen-free F\u00e5hraeus medium and grown in "
        f"controlled-environment chambers (16 h light/8 h dark, 22 C, 55% RH, "
        f"200 umol m-2 s-1); seedlings transferred to sterile 1:1 perlite:vermiculite "
        f"pots and inoculated with Sinorhizobium meliloti 2011 {strain}. Root nodules "
        f"harvested at 3 weeks post-inoculation. mRNA library (oligo-dT), 150 bp "
        f"paired-end, DNBSEQ-T7. Library name: {r['library_name']}. "
        f"Biological replicate {r['replicate']}."
    )
    return title, desc


def write_study_text() -> None:
    txt = f"""# Webin study registration — Medicago DNF2-SRP nodule RNA-seq

- **Study alias:** `{STUDY_ALIAS}`
- **Hold date:** `{HOLD_DATE}` (private until then; release earlier from Webin when the paper is accepted)
- **Submitter account:** Webin-69760
- **Technical contact:** Chongjing Xia <xiachongjing@gmail.com>

## Title

A cross-kingdom interactome predicted by AlphaFold3 reveals a DNF2-centered interface required for symbiotic accommodation

## Abstract

Legumes convert atmospheric nitrogen into ammonium through symbiotic bacteria housed in root nodules, yet the molecular interactions between rhizobial and host proteins inside nodules remain poorly understood. Here we employed AlphaFold3 to construct a cross-kingdom interactome between Medicago truncatula and its symbiont Sinorhizobium meliloti. Screening more than 217,000 protein pairs yielded 7,137 putative interactions, providing a valuable resource for the broader symbiosis community. Within this network, we focused on DEFECTIVE IN NITROGEN FIXATION 2 (DNF2), a host protein required for rhizobial persistence within nodules. We showed that DNF2 localizes to the peribacteroid space and associates with previously uncharacterized secreted rhizobial proteins (SRPs), suggesting it may function as a hub for host-symbiont communication. Notably, knockout of two DNF2-interacting proteins, SRP86 and SRP485, results in white, nitrogen-fixation-deficient nodules with abnormal symbiosomes and elevated expression of senescence-associated genes, closely phenocopying the dnf2 loss-of-function mutant. Together, our findings define a DNF2-SRP molecular framework underlying symbiotic accommodation, and illustrate the potential of AI-guided interactome mapping to uncover molecular mechanisms of plant-microbe interactions with relevance to sustainable agriculture.

## Publication

Preprint: https://www.biorxiv.org/content/10.64898/2026.09.04.749384v1 (DOI 10.64898/2026.09.04.749384)

## Data summary

This ENA study contains 21 paired-end mRNA-seq libraries of Medicago truncatula (ecotype A17) root nodules harvested at 3 weeks post-inoculation: wild type and dnf2 host plants, each inoculated with Sinorhizobium meliloti 2011 wild type or one of five SRP-deletion strains (d-rsp86, d-rsp256, d-rsp485, d-rsp86 d-rsp256, d-rsp86 d-rsp485), three biological replicates per condition. Sequenced 150 bp paired-end on MGI DNBSEQ-T7.
"""
    (ROOT / "study_text.md").write_text(txt)
    print("Wrote study_text.md")


def symlink_fastqs(rows: list[dict]) -> None:
    made = 0
    for r in rows:
        for side in ("fq1", "fq2"):
            fn = r[side]
            dst = ROOT / fn
            if dst.exists() or dst.is_symlink():
                dst.unlink()
            os.symlink(RAW / fn, dst)
            made += 1
    print(f"Symlinked {made} FASTQs into {ROOT.name}/")


def main() -> None:
    md5map = parse_md5file()
    if len(md5map) != 42:
        raise SystemExit(f"expected 42 FASTQ MD5 entries, got {len(md5map)}")
    rows = build_rows(md5map)
    if len(rows) != 21:
        raise SystemExit(f"expected 21 libraries, got {len(rows)}")
    write_metadata(rows)
    write_checklist(rows)
    write_study_text()
    symlink_fastqs(rows)
    print("done")


if __name__ == "__main__":
    main()