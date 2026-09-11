#!/usr/bin/env python3
"""Write experiment/run/submission XML per library (PAIRED, DNBSEQ-T7).

Reads metadata.tsv (which carries the ENA study/sample accessions once filled).
Remote path convention: webin-cli/reads/<LIBRARY_NAME>/<FASTQ>.
"""
from __future__ import annotations

import csv
import html
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "xml" / "runs"

INSTRUMENT = "DNBSEQ-T7"
DRY = os.environ.get("ENA_DRYRUN") == "1"


def esc(s: str) -> str:
    return html.escape(s, quote=True)


def write_one(r: dict) -> None:
    lib = r["library_name"]
    alias = f"webin-reads-{lib}"
    d = OUT / lib
    d.mkdir(parents=True, exist_ok=True)
    study = r.get("study_accession") or ("PRJEB_DRYRUN" if DRY else "")
    sample = r.get("sample_accession") or ("ERS_DRYRUN" if DRY else "")
    remote1 = f"webin-cli/reads/{lib}/{r['fq1']}"
    remote2 = f"webin-cli/reads/{lib}/{r['fq2']}"

    if not study or not sample:
        raise SystemExit(
            f"missing study/sample accession for {lib}; run registration first, "
            f"then fill accessions in metadata.tsv (fill_accessions.py)"
        )

    (d / "experiment.xml").write_text(
        f"""<?xml version="1.0" encoding="UTF-8"?>
<EXPERIMENT_SET>
  <EXPERIMENT alias="{esc(alias)}">
    <TITLE>Raw reads: {esc(lib)}</TITLE>
    <STUDY_REF accession="{esc(study)}"/>
    <DESIGN>
      <DESIGN_DESCRIPTION>150 bp paired-end mRNA-seq, MGI DNBSEQ-T7.</DESIGN_DESCRIPTION>
      <SAMPLE_DESCRIPTOR accession="{esc(sample)}"/>
      <LIBRARY_DESCRIPTOR>
        <LIBRARY_NAME>{esc(lib)}</LIBRARY_NAME>
        <LIBRARY_STRATEGY>RNA-Seq</LIBRARY_STRATEGY>
        <LIBRARY_SOURCE>TRANSCRIPTOMIC</LIBRARY_SOURCE>
        <LIBRARY_SELECTION>Oligo-dT</LIBRARY_SELECTION>
        <LIBRARY_LAYOUT>
          <PAIRED/>
        </LIBRARY_LAYOUT>
      </LIBRARY_DESCRIPTOR>
    </DESIGN>
    <PLATFORM>
      <DNBSEQ>
        <INSTRUMENT_MODEL>{INSTRUMENT}</INSTRUMENT_MODEL>
      </DNBSEQ>
    </PLATFORM>
  </EXPERIMENT>
</EXPERIMENT_SET>
"""
    )
    (d / "run.xml").write_text(
        f"""<?xml version="1.0" encoding="UTF-8"?>
<RUN_SET>
  <RUN alias="{esc(alias)}">
    <TITLE>Raw reads: {esc(lib)}</TITLE>
    <EXPERIMENT_REF refname="{esc(alias)}"/>
    <DATA_BLOCK>
      <FILES>
        <FILE filename="{esc(remote1)}" filetype="fastq" checksum_method="MD5" checksum="{r['fq1_md5']}"/>
        <FILE filename="{esc(remote2)}" filetype="fastq" checksum_method="MD5" checksum="{r['fq2_md5']}"/>
      </FILES>
    </DATA_BLOCK>
  </RUN>
</RUN_SET>
"""
    )
    (d / "submission.xml").write_text(
        f"""<?xml version="1.0" encoding="UTF-8"?>
<SUBMISSION alias="submit-reads-{esc(lib)}">
  <ACTIONS>
    <ACTION>
      <ADD/>
    </ACTION>
  </ACTIONS>
</SUBMISSION>
"""
    )


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    n = 0
    with (ROOT / "metadata.tsv").open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    for r in rows:
        if "TO_COMPUTE" in r["fq1_md5"] or "TO_COMPUTE" in r["fq2_md5"]:
            raise SystemExit(f"missing MD5 for {r['library_name']}")
        write_one(r)
        n += 1
    print(f"Wrote {n} experiment/run/submission triples under {OUT}")


if __name__ == "__main__":
    main()