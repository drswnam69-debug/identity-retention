#!/usr/bin/env python3
"""make_mock_array.py -- synthetic GSE164760-shaped array cohort.

Exists only to exercise the PHASE C code paths (tissue axis, paired
sensitivity, field-effect contrast) before the real data lands. No biology.
"""
from __future__ import annotations

import gzip
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rsi_config as cfg   # noqa: E402

RNG = np.random.default_rng(164760)
GROUPS = [("Healthy liver", 6, 0.0), ("NASH liver", 74, 0.35),
          ("Cirrhotic liver", 8, 0.45),
          ("Non-tumoral NASH liver adjacent to HCC", 29, 0.60),
          ("NASH-HCC tumor", 53, 1.10)]


def main(outdir: str = "data/mock_array") -> int:
    os.makedirs(outdir, exist_ok=True)
    tissue, accs, srcs, eff = [], [], [], []
    n = 0
    tumor_src = 5301
    for label, count, shift in GROUPS:
        for k in range(count):
            tissue.append(label)
            accs.append(f"GSM{5018268 + n}")
            eff.append(shift)
            if label == "NASH-HCC tumor":
                srcs.append(f"NY{tumor_src}")
                tumor_src += 2
            elif label.startswith("Non-tumoral"):
                # first 21 adjacent samples pair with a tumor at src-1
                srcs.append(f"NY{5301 + 2 * k + 1}" if k < 21
                            else f"CF{100 + k}")
            else:
                srcs.append(f"C{n}")
            n += 1

    panel = sorted({g for gs in cfg.ALL_MODULES.values() for g in gs}
                   | set(cfg.SENSITIVITY_S2_CHOLESTEROL))
    genes = panel + [f"FILLER{i:05d}" for i in range(3000)]
    shift = np.array(eff, float)
    base = RNG.uniform(4.0, 11.0, len(genes))[:, None]
    mat = base + RNG.normal(0, 0.7, (len(genes), n))
    for i, g in enumerate(genes):
        if g in cfg.MODULE_SUPPLY or g in cfg.MODULE_REDUCTION:
            mat[i] += 0.8 * shift
        elif g in cfg.MODULE_DRAIN:
            mat[i] -= 0.5 * shift
    # already log-scale, like an RMA-normalized array matrix
    pd.DataFrame(mat, index=genes, columns=accs).to_csv(
        os.path.join(outdir, "mock_array_matrix.tsv.gz"), sep="\t")

    def row(tag, vals):
        return tag + "\t" + "\t".join(f'"{v}"' for v in vals)

    lines = [
        row("!Sample_title", [f"Liver_{t.split()[0]}_{i}"
                              for i, t in enumerate(tissue)]),
        row("!Sample_geo_accession", accs),
        row("!Sample_source_name_ch1", srcs),
        row("!Sample_characteristics_ch1", [f"tissue: {t}" for t in tissue]),
        "!series_matrix_table_begin",
    ]
    with gzip.open(os.path.join(outdir, "mock_array_series.txt.gz"), "wt") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"mock array cohort: {n} samples, {len(genes)} probes -> {outdir}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
