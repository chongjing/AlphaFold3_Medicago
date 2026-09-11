#!/usr/bin/env python3
"""Fill PRJEB/ERS accessions into metadata.tsv from the registration receipts."""
from __future__ import annotations

import csv
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(__file__).resolve().parent
STUDY_REC = ROOT / "xml" / "receipt_prod_study.xml"
SAMPLE_REC = ROOT / "xml" / "receipt_prod_sample.xml"


def study_accession() -> str:
    root = ET.parse(STUDY_REC).getroot()
    for e in root.iter("PROJECT"):
        a = e.get("accession")
        if a:
            return a
    raise SystemExit("no PROJECT accession in study receipt")


def sample_map() -> dict[str, str]:
    root = ET.parse(SAMPLE_REC).getroot()
    m: dict[str, str] = {}
    for s in root.iter("SAMPLE"):
        alias = s.get("alias")
        acc = s.get("accession")
        if alias and acc:
            m[alias] = acc
    return m


def main() -> None:
    study = study_accession()
    sm = sample_map()
    if len(sm) != 21:
        raise SystemExit(f"expected 21 sample accessions, got {len(sm)}")
    path = ROOT / "metadata.tsv"
    rows = list(csv.DictReader(path.open(), delimiter="\t"))
    for row in rows:
        alias = row["sample_alias"]
        if alias not in sm:
            raise SystemExit(f"missing sample accession for alias {alias}")
        row["study_accession"] = study
        row["sample_accession"] = sm[alias]
    fields = list(rows[0].keys())
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"Filled study={study} and {len(sm)} sample accessions into metadata.tsv")


if __name__ == "__main__":
    main()