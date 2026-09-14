#!/usr/bin/env python3
"""20_gse164760_composition_prep.py -- add the 6j composition genes to GSE164760.

The genome-wide GPL13667 probeset matrix is already on disk. Quantile
normalization runs over the COMPLETE probeset matrix before identifier mapping,
so extending the annotation map to carry the composition genes cannot alter the
panel genes' values. That is asserted numerically against the already-written
expr_log, not assumed -- the same check 6h required for the differentiation
panel.
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import io_utils

SRC = "data/GSE164760_probeset.tsv.gz"
PANEL_MAP = "data/GPL13667_combined_map.tsv"
COMP_MAP = "data/GPL13667_composition_probeset_to_symbol.tsv"
COHORT = "results/GSE164760_diff"
OUT = "data/GSE164760_composition_panel.tsv"

print("=== GSE164760 — composition panel (PREREG 6j) ===")
mat = pd.read_csv(SRC, sep="\t", index_col=0)
print(f"  probeset matrix {mat.shape[0]} x {mat.shape[1]}")

# The probeset matrix is keyed by GEO source name; the cohort is keyed by GSM.
# The mapping is the one build_gse164760.py already fixed: SAMPLES in GEO order
# against consecutive accessions from GSM5018268. Read it from that file rather
# than restating it here, so the two cannot drift apart.
_src = open("tools/build_gse164760.py", encoding="utf-8").read()
SAMPLES = _src.split('SAMPLES = """')[1].split('"""')[0].split()
assert len(SAMPLES) == 170, len(SAMPLES)
src2gsm = {s[1:]: f"GSM{5018268 + i}" for i, s in enumerate(SAMPLES)}
assert len(src2gsm) == 170, "duplicate source names"
missing = [c for c in mat.columns if c not in src2gsm]
assert not missing, f"columns with no accession: {missing[:5]}"
mat.columns = [src2gsm[c] for c in mat.columns]

# identical arithmetic to 01_prepare.py --units array-linear
expr = np.log2(mat.clip(lower=0) + 1.0)
expr = io_utils.quantile_normalize(expr)
print("  log2(x+1) then quantile normalization over the complete matrix")


def collapse(mapping_path):
    m = pd.read_csv(mapping_path, sep="\t", header=None,
                    names=["probe", "symbol"])
    have = [p for p in m["probe"] if p in expr.index]
    miss = sorted(set(m["probe"]) - set(have))
    if miss:
        print(f"    !! probesets absent from the matrix: {miss}")
    sub = expr.loc[have].copy()
    sub.index = m.set_index("probe").loc[have, "symbol"].to_numpy()
    return sub.groupby(level=0).max()


panel = collapse(PANEL_MAP)
comp = collapse(COMP_MAP)
print(f"  panel genes {panel.shape[0]}   composition genes {comp.shape[0]}")

# --- assertion: the 6j extension does not move a single panel value ---------
old = pd.read_csv(os.path.join(COHORT, "expr_log.tsv.gz"), sep="\t",
                  index_col=0)
shared = sorted(set(old.index) & set(panel.index))
assert shared, "no overlap with the already-written expr_log"
a = old.loc[shared]
b = panel.loc[shared, old.columns]
d = float(np.abs(a.to_numpy(float) - b.to_numpy(float)).max())
print(f"  [assert] {len(shared)} panel genes reproduce the written expr_log; "
      f"max abs difference {d:.3e}")
assert d < 1e-9, "the composition extension moved panel values"

overlap = sorted(set(comp.index) & set(old.index))
assert not overlap, f"composition genes collide with the panel: {overlap}"

comp = comp[old.columns]
comp.index.name = "gene"
comp.to_csv(OUT, sep="\t")
print(f"  wrote -> {OUT}  ({comp.shape[0]} genes x {comp.shape[1]} samples)")
