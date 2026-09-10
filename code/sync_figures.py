#!/usr/bin/env python3
"""sync_figures.py -- put each rendered figure where the build and the package expect it.

The figure scripts write descriptive names; the manuscript numbers its figures.
The mapping used to live only in whatever shell command last copied the files,
so a regenerated figure could silently fail to reach the built document. It is
written down here and applied for every format at once.
"""
import os
import shutil

HOME = "/home/claude"
MAP = {
    1: "Figure1_concept",
    2: "Figure1_index_spectrum",
    3: "Figure2_tumor_adjacent",
    4: "Figure3_dissociation",
    5: "Figure4_signature_benchmark",
    6: "Figure5_level_transfer",
    7: "Figure6_covariate_comparison",
    8: "Figure7_uncertainty_and_recovery",
    9: "Figure8_third_cohort",
    10: "Figure9_three_tissues",
    11: "Figure11_gao_proteogenomic",
    "S1": "FigureS1_enumeration",
}
SEARCH = [HOME, f"{HOME}/rsi/results/figures"]


def find(stem, ext):
    for d in SEARCH:
        p = os.path.join(d, stem + ext)
        if os.path.exists(p):
            return p
    return None


def main():
    missing = []
    for n, stem in MAP.items():
        png = find(stem, ".png")
        tif = find(stem, ".tiff") or find(stem, ".tif")
        if not png:
            missing.append(f"Figure{n}: no PNG for {stem}")
            continue
        shutil.copy(png, f"{HOME}/build/Figure{n}.png")
        if tif:
            shutil.copy(tif, f"{HOME}/submission_GigaScience/Figure{n}.tif")
        else:
            missing.append(f"Figure{n}: no TIFF for {stem}")
        print(f"  Figure{n:<3} <- {os.path.basename(png)}")
    if missing:
        raise SystemExit("\n".join(["MISSING:"] + missing))
    print("\nall figures synced")


if __name__ == "__main__":
    main()
