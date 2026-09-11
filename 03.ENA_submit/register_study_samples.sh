#!/usr/bin/env bash
# Register the study and 21 samples via the ENA Webin drop-box REST API.
# Credentials from the environment only. Never store the password here.
#
#   export ENA_WEBIN_USER=Webin-69760
#   export ENA_WEBIN_PASSWORD='...'
#   ./register_study_samples.sh test    # wwwdev, discarded in 24 h
#   ./register_study_samples.sh prod    # production PRJEB / ERS
set -euo pipefail

: "${ENA_WEBIN_USER:?set ENA_WEBIN_USER}"
: "${ENA_WEBIN_PASSWORD:?set ENA_WEBIN_PASSWORD}"

MODE="${1:-}"
if [[ "$MODE" != "test" && "$MODE" != "prod" ]]; then
  echo "Usage: $0 test|prod" >&2
  exit 1
fi

HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"
python3 generate_registration_xml.py

if [[ "$MODE" == "test" ]]; then
  URL="https://wwwdev.ebi.ac.uk/ena/submit/drop-box/submit/"
  OUT="$HERE/xml/receipt_test"
  SAMP="xml/sample_test.xml"
else
  URL="https://www.ebi.ac.uk/ena/submit/drop-box/submit/"
  OUT="$HERE/xml/receipt_prod"
  SAMP="xml/sample.xml"
fi

submit() {
  local label="$1"
  shift
  local rec="${OUT}_${label}.xml"
  echo "=== submitting ${label} to ${URL} ==="
  curl -sS -u "${ENA_WEBIN_USER}:${ENA_WEBIN_PASSWORD}" \
    -F "SUBMISSION=@$1" \
    "${@:2}" \
    "$URL" | tee "$rec"
  echo
  if ! grep -q 'success="true"' "$rec"; then
    echo "FAILED: $label  see $rec" >&2
    exit 1
  fi
  echo "OK: $label"
}

if [[ "$MODE" == "test" ]]; then
  submit study xml/submission_study.xml -F "PROJECT=@xml/project.xml"
  submit sample xml/submission_samples.xml -F "SAMPLE=@${SAMP}"
else
  submit study xml/submission_study.xml -F "PROJECT=@xml/project.xml"
  submit sample xml/submission_samples.xml -F "SAMPLE=@${SAMP}"
fi

echo "Receipts in ${OUT}_*.xml"
echo "Next: python3 fill_accessions.py   (writes PRJEB/ERS into metadata.tsv)"