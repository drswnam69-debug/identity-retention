#!/usr/bin/env python3
"""test_pipeline.py -- end-to-end smoke test on synthetic data.

Run this after any code change. It builds a mock cohort, runs both scripts,
and asserts that the index behaves as specified. It proves the code works;
it proves nothing about biology.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import rsi_config as cfg   # noqa: E402
import phenotype as ph     # noqa: E402


def run(cmd: list[str]) -> None:
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise AssertionError(f"FAILED: {' '.join(cmd)}\n{r.stdout}\n{r.stderr}")


def main() -> int:
    tmp = tempfile.mkdtemp(prefix="rsi_test_")
    data, res = os.path.join(tmp, "data"), os.path.join(tmp, "results")

    run([sys.executable, os.path.join(HERE, "make_mock_cohort.py")])
    import make_mock_cohort as mmc
    mmc.main(data)

    run([sys.executable, os.path.join(HERE, "01_prepare.py"),
         "--counts", os.path.join(data, "mock_counts.tsv.gz"),
         "--series", os.path.join(data, "mock_series_matrix.txt.gz"),
         "--cohort", "TEST", "--outdir", res])
    run([sys.executable, os.path.join(HERE, "02_compute_rsi.py"),
         "--cohort", "TEST", "--outdir", res])

    scores = pd.read_csv(os.path.join(res, "TEST", "rsi.tsv"),
                         sep="\t", index_col=0)
    pheno = pd.read_csv(os.path.join(res, "TEST", "phenotype.tsv"),
                        sep="\t", index_col=0)
    qc = json.load(open(os.path.join(res, "TEST", "qc.json")))
    val = json.load(open(os.path.join(res, "TEST",
                                      "internal_validity.json")))

    checks: list[tuple[str, bool]] = []

    checks.append(("all samples scored",
                   scores.shape[0] == qc["n_samples"] == 216))
    checks.append(("RSI is finite everywhere",
                   bool(np.isfinite(scores["RSI"]).all())))

    # RSI must equal the locked formula exactly
    recomputed = (cfg.WEIGHT_SUPPLY * scores["z_supply"]
                  + cfg.WEIGHT_REDUCTION * scores["z_reduction"]
                  - cfg.WEIGHT_DRAIN * scores["z_drain"])
    checks.append(("RSI matches the locked formula",
                   bool(np.allclose(recomputed, scores["RSI"]))))

    # every module z-score must be mean 0 across the cohort (within-cohort z)
    checks.append(("component z-scores are cohort-centered",
                   bool(np.allclose(scores[["z_supply", "z_reduction",
                                            "z_drain"]].mean(), 0,
                                    atol=1e-9))))

    checks.append(("no specified gene silently missing",
                   all(not v for v in val["missing_genes"].values())))

    expr = pd.read_csv(os.path.join(res, "TEST", "expr_log.tsv.gz"),
                       sep="\t", index_col=0)
    checks.append(("aligner summary rows are excluded from the matrix",
                   not expr.index.astype(str).str.startswith(("__", "N_")).any()))

    checks.append(("stage column fully assigned",
                   qc["n_stage_unassigned"] == 0))
    checks.append(("stage order preserved",
                   list(qc["stage_counts"]) == ph.STAGE_ORDER))

    # the mock data was built with an upward RSI trend; the pipeline must see it
    order = {s: i for i, s in enumerate(ph.STAGE_ORDER)}
    sev = pheno.loc[scores.index, "stage"].map(order)
    # Spearman via numpy ranks, so the smoke test needs no SciPy
    ra = pd.Series(scores["RSI"].values).rank().to_numpy()
    rb = pd.Series(sev.values, dtype=float).rank().to_numpy()
    rho = float(np.corrcoef(ra, rb)[0, 1])
    checks.append((f"recovers the planted trend (rho={rho:.2f} > 0.5)",
                   rho > 0.5))

    checks.append(("lock hash is stable",
                   val["lock_hash"] == cfg.lock_hash()))

    # ---- PHASE B integration, including a negative control -------------
    run([sys.executable, os.path.join(HERE, "03_phase_b.py"),
         "--cohort", "TEST", "--outdir", res])
    pb = json.load(open(os.path.join(res, "TEST", "phase_b.json")))
    pr = pb["primary_RSI"]
    checks.append((f"Phase B detects the planted trend (p={pr['p_increasing']:.2g})",
                   pr["p_increasing"] < 1e-6))
    checks.append(("Phase B stage medians are monotone increasing",
                   all(a < b for a, b in zip(pr["median_per_stage"],
                                             pr["median_per_stage"][1:]))))
    checks.append(("Phase B recovers the planted component directions",
                   pb["components"]["z_supply"]["z"] > 0
                   and pb["components"]["z_reduction"]["z"] > 0
                   and pb["components"]["z_drain"]["z"] < 0))

    # negative control: shuffling the stage labels must destroy the signal
    shuf = os.path.join(res, "SHUFFLE")
    os.makedirs(shuf, exist_ok=True)
    import shutil
    for f in ("rsi.tsv", "expr_log.tsv.gz", "qc.json"):
        shutil.copy(os.path.join(res, "TEST", f), os.path.join(shuf, f))
    ph_df = pd.read_csv(os.path.join(res, "TEST", "phenotype.tsv"), sep="\t",
                        index_col=0)
    ph_df["stage"] = (np.random.default_rng(3)
                      .permutation(ph_df["stage"].to_numpy()))
    ph_df.to_csv(os.path.join(shuf, "phenotype.tsv"), sep="\t")
    run([sys.executable, os.path.join(HERE, "03_phase_b.py"),
         "--cohort", "SHUFFLE", "--outdir", res])
    pn = json.load(open(os.path.join(res, "SHUFFLE", "phase_b.json")))
    p_null = pn["primary_RSI"]["p_increasing"]
    checks.append((f"shuffled stages give no trend (p={p_null:.2g})",
                   p_null > 0.05))

    # ---- PHASE C on a synthetic array cohort ---------------------------
    import make_mock_array as mma
    adir = os.path.join(tmp, "data_array")
    mma.main(adir)
    run([sys.executable, os.path.join(HERE, "01_prepare.py"),
         "--counts", os.path.join(adir, "mock_array_matrix.tsv.gz"),
         "--series", os.path.join(adir, "mock_array_series.txt.gz"),
         "--cohort", "ARRAY", "--outdir", res, "--array",
         "--stage-source", "tissue"])
    run([sys.executable, os.path.join(HERE, "02_compute_rsi.py"),
         "--cohort", "ARRAY", "--outdir", res])
    run([sys.executable, os.path.join(HERE, "04_phase_c.py"),
         "--cohort", "ARRAY", "--outdir", res])
    pc = json.load(open(os.path.join(res, "ARRAY", "phase_c.json")))

    checks.append(("Phase C reads the five GSE164760 tissue groups",
                   list(pc["group_sizes"].values()) == [6, 74, 8, 29, 53]))
    checks.append((f"Phase C detects the planted H2 shift "
                   f"(p={pc['H2']['p_two_sided']:.1g})",
                   pc["H2"]["p_two_sided"] < 0.01
                   and pc["H2"]["hodges_lehmann_shift"] > 0))
    checks.append((f"Phase C detects the planted H3 field effect "
                   f"(p={pc['H3']['p_two_sided']:.1g})",
                   pc["H3"]["p_two_sided"] < 0.01
                   and pc["H3"]["hodges_lehmann_shift"] > 0))
    checks.append((f"paired sensitivity finds the planted pairs "
                   f"(n={pc['H2_paired_sensitivity'].get('n_pairs')})",
                   pc["H2_paired_sensitivity"].get("n_pairs") == 21))
    checks.append(("Cliff's delta stays within [-1, 1]",
                   all(-1 <= pc[k]["cliffs_delta"] <= 1 for k in ("H2", "H3"))))

    # negative control: shuffling tissue labels must kill the field effect
    ashuf = os.path.join(res, "ARRAY_SHUF")
    os.makedirs(ashuf, exist_ok=True)
    for f in ("rsi.tsv", "qc.json"):
        shutil.copy(os.path.join(res, "ARRAY", f), os.path.join(ashuf, f))
    ap = pd.read_csv(os.path.join(res, "ARRAY", "phenotype.tsv"), sep="\t",
                     index_col=0)
    ap["tissue"] = np.random.default_rng(11).permutation(
        ap["tissue"].to_numpy())
    ap.to_csv(os.path.join(ashuf, "phenotype.tsv"), sep="\t")
    run([sys.executable, os.path.join(HERE, "04_phase_c.py"),
         "--cohort", "ARRAY_SHUF", "--outdir", res])
    pcn = json.load(open(os.path.join(res, "ARRAY_SHUF", "phase_c.json")))
    checks.append((f"shuffled tissue labels give no field effect "
                   f"(p={pcn['H3']['p_two_sided']:.2g})",
                   pcn["H3"]["p_two_sided"] > 0.05))

    width = max(len(n) for n, _ in checks)
    ok = True
    for name, passed in checks:
        print(f"  [{'PASS' if passed else 'FAIL'}] {name:<{width}}")
        ok &= passed
    print(f"\n{'ALL CHECKS PASSED' if ok else 'SOME CHECKS FAILED'}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
