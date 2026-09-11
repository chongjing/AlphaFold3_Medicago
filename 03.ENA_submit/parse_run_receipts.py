#!/usr/bin/env python3
"""Merge xml/run_receipts.tsv + metadata.tsv -> accessions_registered.tsv + DAS."""
from __future__ import annotations

import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parent
HOLD = "2026-11-01"
DOI = "10.64898/2026.09.04.749384"
DOI_URL = "https://www.biorxiv.org/content/10.64898/2026.09.04.749384v1"


def main() -> None:
    meta = {r["library_name"]: r
            for r in csv.DictReader((ROOT / "metadata.tsv").open(), delimiter="\t")}
    recs = {}
    with (ROOT / "xml" / "run_receipts.tsv").open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            recs[r["library_name"]] = r

    study = meta[list(meta)[0]]["study_accession"] or "PRJEB_UNKNOWN"

    # accessions_registered.tsv
    out = ROOT / "accessions_registered.tsv"
    cols = ["library_name", "sample_alias", "host_genotype", "strain", "replicate",
            "ERS", "ERR", "ERX", "PRJEB"]
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for lib in meta:  # preserve metadata order
            r = recs.get(lib, {})
            w.writerow({
                "library_name": lib,
                "sample_alias": meta[lib]["sample_alias"],
                "host_genotype": meta[lib]["host_genotype"],
                "strain": meta[lib]["strain"],
                "replicate": meta[lib]["replicate"],
                "ERS": meta[lib].get("sample_accession", ""),
                "ERR": r.get("ERR", ""),
                "ERX": r.get("ERX", ""),
                "PRJEB": study,
            })

    # data_availability_statement.md
    das = ROOT / "data_availability_statement.md"
    lines = [
        "# Data Availability Statement",
        "",
        f"All raw RNA-seq reads generated in this study have been deposited in the",
        f"European Nucleotide Archive (ENA) at EMBL-EBI under project accession",
        f"**{study}** (https://www.ebi.ac.uk/ena/browser/view/{study}).",
        f"Data are held private until **{HOLD}** and can be released earlier from the",
        f"Webin portal when the paper is accepted.",
        "",
        f"This study is associated with the preprint: {DOI_URL} (DOI: {DOI}).",
        "",
        "| Library | Host genotype | Rhizobial strain | Rep | Sample (ERS) | Run (ERR) | Experiment (ERX) |",
        "|---|---|---|---|---|---|---|",
    ]
    for lib in meta:
        m = meta[lib]
        r = recs.get(lib, {})
        lines.append(
            f"| {lib} | {m['host_genotype']} | {m['strain']} | {m['replicate']} "
            f"| {m.get('sample_accession','')} | {r.get('ERR','')} | {r.get('ERX','')} |"
        )
    lines += [
        "",
        "**Accession summary:**",
        f"- BioProject: {study}",
        "- 21 samples, 42 paired-end FASTQ files, ~49 GB",
        "- Instrument: MGI DNBSEQ-T7, 150 bp paired-end mRNA (oligo-dT), BGI",
        f"- Preprint: {DOI_URL} (DOI {DOI})",
        "",
    ]
    das.write_text("\n".join(lines))
    ok = sum(1 for r in recs.values() if r.get("success") == "true")
    print(f"Wrote {out}, {das}; runs ok={ok}/{len(recs)}")


if __name__ == "__main__":
    main()