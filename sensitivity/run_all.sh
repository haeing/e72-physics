#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p results
for stage in yield acceptance diff_cross polar_th crystalball export_pdf; do
  echo "Running ${stage}"
  if [[ "$stage" == export_pdf ]]; then
    root -l -b -q macro/export_pdf.C > "results/${stage}.log" 2>&1
  else
    root -l -b -q "macro/run_stage.C(\"${stage}\")" > "results/${stage}.log" 2>&1
  fi
  if grep -Eq '(^|[[:space:]])(error:|Error in <|fatal error:)|segmentation violation' "results/${stage}.log"; then
    echo "ROOT reported an error; inspect results/${stage}.log" >&2
    exit 1
  fi
done
echo 'Done: results/sensitivity.pdf, results/*.root, results/*.log'
