#!/usr/bin/env python3
"""make_mock_cohort.py -- synthetic GEO-shaped data to test the pipeline.

This exists ONLY to prove the code runs end to end before the real files
arrive. It contains no biology and must never be used for any result.
"""

from __future__ import annotations

import gzip
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rsi_config as cfg   # noqa: E402
import phenotype as ph     # noqa: E402

RNG = np.random.default_rng(20260825)
N_PER_STAGE = {"control": 10, "NAFL": 51, "NASH_F0-F1": 34,
               "NASH_F2": 53, "NASH_F3": 54, "NASH_F4": 14}


def main(outdir: str = "data/mock") -> int:
    os.makedirs(outdir, exist_ok=True)
    stages, samples = [], []
    for stage, n in N_PER_STAGE.items():
        for _ in range(n):
            stages.append(stage)
            samples.append(f"GSM{9000000 + len(samples)}")

    panel = sorted({g for gs in cfg.ALL_MODULES.values() for g in gs}
                   | set(cfg.SENSITIVITY_S2_CHOLESTEROL)
                   | {"XIST", "RPS4Y1", "DDX3Y"})
    filler = [f"FILLER{i:05d}" for i in range(4000)]
    genes = panel + filler

    sev = np.array([ph.STAGE_ORDER.index(s) for s in stages], float)
    sev = (sev - sev.mean()) / sev.std()

    base = RNG.uniform(4.0, 9.0, size=len(genes))[:, None]
    mat = base + RNG.normal(0, 0.6, size=(len(genes), len(samples)))
    for i, g in enumerate(genes):
        if g in cfg.MODULE_SUPPLY or g in cfg.MODULE_REDUCTION:
            mat[i] += 0.45 * sev
        elif g in cfg.MODULE_DRAIN:
            mat[i] -= 0.30 * sev
    counts = pd.DataFrame(
        np.round(2 ** mat).astype(int), index=genes, columns=samples)
    # real GEO count files carry aligner summary rows; they must be dropped
    # before normalization, so put some in the fixture and check for them.
    summary = pd.DataFrame(
        RNG.integers(5e5, 5e6, size=(5, len(samples))),
        index=["__no_feature", "__ambiguous", "__too_low_aQual",
               "__not_aligned", "__alignment_not_unique"],
        columns=samples)
    counts = pd.concat([counts, summary])
    counts.to_csv(os.path.join(outdir, "mock_counts.tsv.gz"), sep="\t")

    lines = ['!Sample_title\t' + "\t".join(f'"{s}_liver"' for s in samples),
             '!Sample_geo_accession\t' + "\t".join(f'"{s}"' for s in samples),
             '!Sample_characteristics_ch1\t'
             + "\t".join(f'"group in paper: {s}"' for s in stages),
             '!Sample_characteristics_ch1\t'
             + "\t".join(f'"nas score: {RNG.integers(0, 9)}"' for _ in stages),
             '!series_matrix_table_begin']
    with gzip.open(os.path.join(outdir, "mock_series_matrix.txt.gz"),
                   "wt") as fh:
        fh.write("\n".join(lines) + "\n")

    print(f"mock cohort: {len(samples)} samples, {len(genes)} genes "
          f"-> {outdir}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
