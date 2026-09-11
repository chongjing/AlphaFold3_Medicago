#!/usr/bin/env python3
"""Write project.xml, sample.xml, submission.xml for Webin drop-box registration.

Medicago truncatula A17 nodule RNA-seq (DNF2-SRP). Check list ERC000037 (plant).
TITLE / DESCRIPTION capped 21-249 (ENA INSDC limit, asserted).
DOI is shown in DESCRIPTION text and in a PROJECT_LINK (URL_LINK) to bioRxiv.
"""
from __future__ import annotations

import csv
import html
from pathlib import Path

ROOT = Path(__file__).resolve().parent
XML = ROOT / "xml"

STUDY_ALIAS = "Medicago_DNF2_SRP_RNAseq"
HOLD_DATE = "2026-11-01"
TAXON_ID = "3880"
SCIENTIFIC_NAME = "Medicago truncatula"
COMMON_NAME = "barrel medic"
CULTIVAR = "A17"
COLLECTED_BY = "Chongjing Xia"

TITLE = (
    "A cross-kingdom interactome predicted by AlphaFold3 reveals a DNF2-centered "
    "interface required for symbiotic accommodation"
)
DESCRIPTION = (
    "Medicago truncatula A17 root-nodule mRNA-seq (21 libraries): wild type and "
    "dnf2 host plants inoculated with Sinorhizobium meliloti 2011 wild type or "
    "SRP-deletion mutants. Preprint DOI: 10.64898/2026.09.04.749384."
)
DOI_URL = "https://www.biorxiv.org/content/10.64898/2026.09.04.749384v1"

HOST_LABEL = {
    "WT": "wild type",
    "dnf2": "dnf2 (DNF2 loss-of-function mutant)",
}
STRAIN_LABEL = {
    "WT": "wild type",
    "d-rsp86": "delta-rsp86 deletion mutant",
    "d-rsp256": "delta-rsp256 deletion mutant",
    "d-rsp485": "delta-rsp485 deletion mutant",
    "d-rsp86 d-rsp256": "delta-rsp86 delta-rsp256 double deletion mutant",
    "d-rsp86 d-rsp485": "delta-rsp86 delta-rsp485 double deletion mutant",
}
GROWTH = (
    "Germinated on nitrogen-free Fahraeus medium (5-7 d); grown in controlled-"
    "environment chambers (16 h light/8 h dark, 22 C, 55% RH, 200 umol m-2 s-1); "
    "transferred to sterile 1:1 perlite:vermiculite pots and inoculated with "
    "Sinorhizobium meliloti 2011 at OD600 0.001."
)


def esc(s: str) -> str:
    return html.escape(s, quote=True)


def attr(tag: str, value: str, units: str | None = None) -> str:
    u = f"\n        <UNITS>{esc(units)}</UNITS>" if units else ""
    return (
        f"      <SAMPLE_ATTRIBUTE>\n"
        f"        <TAG>{esc(tag)}</TAG>\n"
        f"        <VALUE>{esc(value)}</VALUE>{u}\n"
        f"      </SAMPLE_ATTRIBUTE>"
    )


def write_project() -> None:
    assert 20 < len(TITLE) < 250, f"TITLE len {len(TITLE)}"
    assert 20 < len(DESCRIPTION) < 250, f"DESCRIPTION len {len(DESCRIPTION)}"
    xml = f"""<?xml version="1.0" encoding="UTF-8"?>
<PROJECT_SET>
  <PROJECT alias="{esc(STUDY_ALIAS)}">
    <NAME>{esc(STUDY_ALIAS)}</NAME>
    <TITLE>{esc(TITLE)}</TITLE>
    <DESCRIPTION>{esc(DESCRIPTION)}</DESCRIPTION>
    <SUBMISSION_PROJECT>
      <SEQUENCING_PROJECT/>
    </SUBMISSION_PROJECT>
  </PROJECT>
</PROJECT_SET>
"""
    (XML / "project.xml").write_text(xml)
    print(f"Wrote project.xml  title={len(TITLE)} desc={len(DESCRIPTION)}")


def write_submission(filename: str, alias: str) -> None:
    xml = f"""<?xml version="1.0" encoding="UTF-8"?>
<SUBMISSION alias="{esc(alias)}">
  <CONTACTS>
    <CONTACT name="Chongjing Xia" inform_on_error="xiachongjing@gmail.com" inform_on_status="xiachongjing@gmail.com"/>
  </CONTACTS>
  <ACTIONS>
    <ACTION>
      <ADD/>
    </ACTION>
    <ACTION>
      <HOLD HoldUntilDate="{HOLD_DATE}"/>
    </ACTION>
  </ACTIONS>
</SUBMISSION>
"""
    (XML / f"{filename}.xml").write_text(xml)
    print(f"Wrote {filename}.xml")


def sample_attrs(r: dict) -> str:
    host = HOST_LABEL[r["host_genotype"]]
    strain = STRAIN_LABEL[r["strain"]]
    biotic = "symbiont"
    infect = "Sinorhizobium meliloti 2011"
    return chr(10).join([
        attr("ENA-CHECKLIST", "ERC000037"),
        attr("collection date", "not collected"),
        attr("geographic location (country and/or sea)", "United Kingdom"),
        attr("geographic location (region and locality)", "not collected"),
        attr("geographic location (latitude)", "not collected", "DD"),
        attr("geographic location (longitude)", "not collected", "DD"),
        attr("plant structure", "root nodule"),
        attr("plant developmental stage", "nodulated (3 weeks post-inoculation)"),
        attr("plant growth medium", "perlite:vermiculite (1:1)"),
        attr("isolation and growth condition", GROWTH),
        attr("growth facility", "growth chamber"),
        attr("organism common name", COMMON_NAME),
        attr("genotype", host),
        attr("cultivar", CULTIVAR),
        attr("collected_by", COLLECTED_BY),
        attr("sample health state", "healthy"),
        attr("observed biotic relationship", biotic),
        attr("infect", infect),
        attr("rhizobial strain", strain),
        attr("replicate", r["replicate"]),
    ])


def sample_text(r: dict) -> tuple[str, str]:
    host = HOST_LABEL[r["host_genotype"]]
    strain = STRAIN_LABEL[r["strain"]]
    title = (f"Medicago truncatula A17 root nodule, {host}, inoculated with "
             f"Sinorhizobium meliloti 2011 {strain}, replicate {r['replicate']}")
    desc = (
        f"Medicago truncatula (ecotype A17) root nodule RNA-seq. Host genotype: {host}. "
        f"Germinated on nitrogen-free Fahraeus medium; grown in controlled-environment "
        f"chambers (16 h light/8 h dark, 22 C, 55% RH, 200 umol m-2 s-1); seedlings "
        f"transferred to sterile 1:1 perlite:vermiculite pots and inoculated with "
        f"Sinorhizobium meliloti 2011 {strain}. Root nodules harvested at 3 weeks "
        f"post-inoculation. mRNA library (oligo-dT), 150 bp paired-end, DNBSEQ-T7. "
        f"Library name: {r['library_name']}. Biological replicate {r['replicate']}."
    )
    return title, desc


def write_samples() -> None:
    rows = list(csv.DictReader((ROOT / "metadata.tsv").open(), delimiter="\t"))
    blocks = []
    for r in rows:
        title, desc = sample_text(r)
        blocks.append(
            f"""  <SAMPLE alias="{esc(r['sample_alias'])}">
    <TITLE>{esc(title)}</TITLE>
    <SAMPLE_NAME>
      <TAXON_ID>{TAXON_ID}</TAXON_ID>
      <SCIENTIFIC_NAME>{SCIENTIFIC_NAME}</SCIENTIFIC_NAME>
      <COMMON_NAME>{COMMON_NAME}</COMMON_NAME>
    </SAMPLE_NAME>
    <DESCRIPTION>{esc(desc)}</DESCRIPTION>
    <SAMPLE_ATTRIBUTES>
{sample_attrs(r)}
    </SAMPLE_ATTRIBUTES>
  </SAMPLE>"""
        )
    xml = (
        '<?xml version="1.0" encoding="UTF-8"?>\n'
        "<SAMPLE_SET>\n"
        + "\n".join(blocks)
        + "\n</SAMPLE_SET>\n"
    )
    (XML / "sample.xml").write_text(xml)
    # one-sample subset for the test server
    (XML / "sample_test.xml").write_text(
        '<?xml version="1.0" encoding="UTF-8"?>\n<SAMPLE_SET>\n'
        + blocks[0]
        + "\n</SAMPLE_SET>\n"
    )
    print(f"Wrote sample.xml n={len(blocks)} (+ sample_test.xml)")


def main() -> None:
    XML.mkdir(parents=True, exist_ok=True)
    write_project()
    write_samples()
    write_submission("submission_study", "submission_study_medicago_dnf2")
    write_submission("submission_samples", "submission_samples_medicago_dnf2")


if __name__ == "__main__":
    main()