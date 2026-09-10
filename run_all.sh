#!/usr/bin/env bash
# run_all.sh -- PHASE A and PHASE B for every cohort present in data/
set -euo pipefail
cd "$(dirname "$0")"

echo "=== verifying the statistics ==="
python3 code/test_stats.py
echo
echo "=== verifying the pipeline on synthetic data ==="
python3 code/test_pipeline.py

DONE=()

run_cohort () {  # name  counts  series  [extra flags...]
  local name="$1" counts="$2" series="$3"; shift 3
  if [[ ! -f "$counts" || ! -f "$series" ]]; then
    echo "--- skip $name (files not in data/) ---"; return
  fi
  echo; echo "=== $name ==="
  python3 code/01_prepare.py --counts "$counts" --series "$series" \
      --cohort "$name" --outdir results "$@"
  python3 code/02_compute_rsi.py --cohort "$name" --outdir results
  python3 code/02_compute_rsi.py --cohort "$name" --outdir results \
      --sensitivity-cholesterol
  python3 code/03_phase_b.py --cohort "$name" --outdir results
  DONE+=("$name")
}

ANNOT_ARGS=()
[[ -f data/Human.GRCh38.p13.annot.tsv.gz ]] && \
  ANNOT_ARGS=(--annotation data/Human.GRCh38.p13.annot.tsv.gz)

run_cohort GSE135251 data/GSE135251_RAW.tar \
                     data/GSE135251_series_matrix.txt.gz "${ANNOT_ARGS[@]}"
run_cohort GSE130970 data/GSE130970_all_sample_salmon_tximport_counts_entrez_gene_ID.csv.gz \
                     data/GSE130970_series_matrix.txt.gz "${ANNOT_ARGS[@]}"
run_cohort GSE167523 data/GSE167523_counts.txt.gz \
                     data/GSE167523_series_matrix.txt.gz
run_cohort GSE126848 data/GSE126848_Gene_counts_raw.txt.gz \
                     data/GSE126848_series_matrix.txt.gz

if [[ ${#DONE[@]} -gt 0 ]]; then
  echo; echo "=== Figure 2 ==="
  python3 code/05_figures.py --cohorts "${DONE[@]}" \
      --outdir results --figdir results/figures
else
  echo; echo "No cohorts found. See DATA_ACQUISITION.md."
fi

echo; echo "PHASE A and B complete. Outcomes remain unopened."
