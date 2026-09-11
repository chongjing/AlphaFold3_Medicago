# 03.ENA_submit — ENA raw-read deposition (PRJEB126112)

Submission of the **Medicago truncatula DNF2–SRP nodule RNA-seq** (batch J003) to the
European Nucleotide Archive (ENA) at EMBL-EBI, via the Webin drop-box REST API + curl FTP.

## Accessions (PRIVATE until 2026-11-01)

| Object | Accession |
|---|---|
| BioProject | [PRJEB126112](https://www.ebi.ac.uk/ena/browser/view/PRJEB126112) (study ERP205777) |
| Samples | ERS31306659 – ERS31306679 (21) |
| Experiments | ERX17259589 – ERX17259609 (21) |
| Runs | ERR17869048 – ERR17869068 (21) |

Full mapping: [`accessions_registered.tsv`](accessions_registered.tsv)
Manuscript data-availability text: [`data_availability_statement.md`](data_availability_statement.md)

## What was deposited

- 21 paired-end mRNA-seq libraries — 7 conditions × 3 biological replicates
- _Medicago truncatula_ (ecotype A17) root nodules, 3 weeks post-inoculation
- Host genotypes: WT and `dnf2`; rhizobial strains: Sm2011 WT and five SRP-deletion mutants
  (`Δrsp86`, `Δrsp256`, `Δrsp485`, `Δrsp86 Δrsp256`, `Δrsp86 Δrsp485`)
- Instrument: MGI **DNBSEQ-T7**, 150 bp paired-end, polyA (oligo-dT)
- 42 FASTQ files, ~49 GB (data live on ENA, not in this repo)
- Preprint: https://www.biorxiv.org/content/10.64898/2026.09.04.749384v1
  (DOI `10.64898/2026.09.04.749384`)

## Files

| File | Role |
|---|---|
| `generate_ena_pack.py` | build metadata.tsv + ERC000037 checklist + study text + FASTQ symlinks |
| `generate_registration_xml.py` | study + sample XML (drop-box REST) |
| `generate_run_xml.py` | experiment/run/submission XML per library (PAIRED, DNBSEQ) |
| `fill_accessions.py` | write PRJEB/ERS back into metadata.tsv after registration |
| `upload_fastq.py` / `upload_fastq.srun` | resumable curl-FTP upload (state-file) |
| `register_study_samples.sh` | register study + 21 samples (test/prod) |
| `submit_runs_rest.sh` | register experiments + runs |
| `parse_run_receipts.py` | merge receipts into accessions + data-availability statement |
| `metadata.tsv` | master table (21 libraries, MD5, byte sizes) |
| `sample_checklist_ERC000037.tsv` | plant checklist (ERC000037) spreadsheet |
| `study_text.md` | study title / abstract / DOI / contact |
| `accessions_registered.tsv` | final sample → ERS/ERR/ERX mapping |
| `data_availability_statement.md` | manuscript-ready data-availability paragraph |

## Reproduce / notes

- Credentials are passed via environment (`ENA_WEBIN_USER`, `ENA_WEBIN_PASSWORD`) — never stored here.
- The FASTQ files are on ENA; this directory holds the reproducible submission scripts + metadata
  only (symlinks to raw data are intentionally excluded from the repo).
- Platform block in `experiment.xml` is `<DNBSEQ><INSTRUMENT_MODEL>DNBSEQ-T7</INSTRUMENT_MODEL></DNBSEQ>`
  with an explicit `<PAIRED/>` layout (webin-cli otherwise mis-writes `<SINGLE/>`).