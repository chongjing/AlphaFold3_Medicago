#!/bin/bash
# ===========================================================================
# Run the unified edgeR DE analysis (all genotypes vs WT).
#   Usage: nohup ./01.run_DE.sh 1>01.run_DE.log 2>&1 &
#   Subset: ./01.run_DE.sh A256 dnf2   # only these genotypes
# ===========================================================================
set -euo pipefail

RSCRIPT="Rscript"
SCRIPT_DIR="../05.DE_analysis"
SCRIPT="${SCRIPT_DIR}/00.edgeR_unified.R"

echo "=== edgeR unified DE analysis vs WT ==="
echo "Started: $(date)"
cd "$SCRIPT_DIR"

"$RSCRIPT" "$SCRIPT" "$@"

echo "Finished: $(date)"
echo ""
echo "=== Summary ==="
sum="${SCRIPT_DIR}/99.Summary/All_comparisons_DE_summary.tsv"
[ -f "$sum" ] && column -t -s$'\t' "$sum"
