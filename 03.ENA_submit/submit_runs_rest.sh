#!/usr/bin/env bash
# After FASTQs are on webin2, register experiment+run via drop-box REST (PAIRED, DNBSEQ).
set -euo pipefail

: "${ENA_WEBIN_USER:?set ENA_WEBIN_USER}"
: "${ENA_WEBIN_PASSWORD:?set ENA_WEBIN_PASSWORD}"

HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"
python3 generate_run_xml.py
URL="https://www.ebi.ac.uk/ena/submit/drop-box/submit/"
SUMMARY="$HERE/xml/run_receipts.tsv"
echo -e "library_name\tsuccess\tERR\tERX\tERA\treceipt" > "$SUMMARY"

fail=0
for d in "$HERE"/xml/runs/*/; do
  lib="$(basename "$d")"
  rec="$d/receipt.xml"
  echo "=== REST submit $lib ==="
  curl -sS -u "${ENA_WEBIN_USER}:${ENA_WEBIN_PASSWORD}" \
    -F "SUBMISSION=@${d}submission.xml" \
    -F "EXPERIMENT=@${d}experiment.xml" \
    -F "RUN=@${d}run.xml" \
    "$URL" | tee "$rec"
  echo
  if grep -q 'success="true"' "$rec"; then
    acc=$(python3 - "$rec" <<'PY'
import sys, xml.etree.ElementTree as ET
r = ET.parse(sys.argv[1]).getroot()
run = r.find("RUN"); exp = r.find("EXPERIMENT"); sub = r.find("SUBMISSION")
print((run.get("accession") if run is not None else ""),
      (exp.get("accession") if exp is not None else ""),
      (sub.get("accession") if sub is not None else ""))
PY
)
    set -- $acc
    echo -e "${lib}\ttrue\t${1:-}\t${2:-}\t${3:-}\t${rec}" >> "$SUMMARY"
    echo "OK $lib ERR=${1:-} ERX=${2:-}"
  else
    echo -e "${lib}\tfalse\t\t\t\t${rec}" >> "$SUMMARY"
    echo "FAILED $lib  see $rec" >&2
    fail=1
  fi
done

echo "Summary: $SUMMARY"
exit $fail