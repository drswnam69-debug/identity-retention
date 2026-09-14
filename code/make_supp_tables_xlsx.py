#!/usr/bin/env python3
"""make_supp_tables_xlsx.py -- build Supplementary Material 2.

Supplementary Tables S1 to S9 used to reach the editor as twelve CSV files
inside the code deposit, which the submission system unpacked into the several
hundred files of that deposit. The code deposit is now cited by its Zenodo
digital object identifier instead, and the tables, which are the part of it a
reader of the article actually needs, are gathered here into a single workbook
with one worksheet per file and a contents sheet in front.

Nothing is recomputed. Every cell is copied from the archived CSV, and the
contents sheet records which file each worksheet came from, so the workbook can
be rebuilt from the deposit and checked against it.
"""
import csv
import os

from openpyxl import Workbook
from openpyxl.styles import Alignment, Font
from openpyxl.utils import get_column_letter

SRC = "/home/claude/SF2build/SupplementaryFile2"
OUT = "/home/claude/build/04_Supplementary Material 2_Supplementary tables.xlsx"

# worksheet name -> (source file, one-line description)
SHEETS = [
    ("S1_signature_retention", "SupplementaryTable_S1_signature_retention.csv",
     "All 126 size-eligible liver signatures: gene counts, covariate overlaps, "
     "unadjusted shifts, retention before and after removing shared genes."),
    ("S2_GSE14520", "SupplementaryTable_S2_covariate_comparison_GSE14520.csv",
     "GSE14520: retention under identity, under composition and under the joint "
     "model, with the signed joint intercept and the flags."),
    ("S2_GSE76427", "SupplementaryTable_S2_covariate_comparison_GSE76427.csv",
     "GSE76427: the same quantities as S2_GSE14520."),
    ("S3_GSE14520", "SupplementaryTable_S3_retention_intervals_GSE14520.csv",
     "GSE14520: 95% paired-bootstrap interval per signature, interval width and "
     "instability flags."),
    ("S3_GSE76427", "SupplementaryTable_S3_retention_intervals_GSE76427.csv",
     "GSE76427: the same quantities as S3_GSE14520."),
    ("S4_simulation_grid", "SupplementaryTable_S4_simulation_grid.csv",
     "All 34 pre-registered simulation cells: median estimate, bias, "
     "interquartile range and achieved coverage."),
    ("S5_TCGA_LIHC", "SupplementaryTable_S5_TCGA_LIHC_retention.csv",
     "The 103 evaluable TCGA-LIHC signatures, with the signed intercept and flags."),
    ("S6_kidney", "SupplementaryTable_S6_kidney_retention.csv",
     "The evaluable TCGA-KIRC signatures, with the signed intercept and flags."),
    ("S6_lung", "SupplementaryTable_S6_lung_retention.csv",
     "The evaluable TCGA-LUAD signatures, with the signed intercept and flags."),
    ("S7_prognostic", "SupplementaryTable_S7_prognostic.csv",
     "Every evaluable TCGA-LIHC signature: hazard ratio before and after "
     "adjustment for the tumor's identity score, Benjamini-Hochberg q, and the "
     "surviving fraction of the log hazard ratio."),
    ("S8_enumeration", "SupplementaryTable_S8_enumeration.csv",
     "Every enumerated gene set in the three tissues, 223 rows: size, symbols "
     "present on the platform, unadjusted paired shift and P value, whether it "
     "entered the retention distribution, and the reason for any exclusion."),
    ("S9_negative_control", "SupplementaryTable_S9_negative_control.csv",
     "The 119 GSE14520 signatures: retention under D1, under the quadratic model "
     "of 6ae, and the median and 95% range across the 200 matched random "
     "covariates."),
]

BODY = Font(name="Arial", size=10)
HEAD = Font(name="Arial", size=10, bold=True)
TITLE = Font(name="Arial", size=12, bold=True)


def num(v):
    """Numbers go in as numbers so the workbook can be sorted and plotted."""
    if v == "":
        return None
    try:
        return int(v) if v.lstrip("-").isdigit() else float(v)
    except ValueError:
        return v


def write_sheet(ws, rows):
    for r in rows:
        ws.append([num(c) for c in r])
    for c in ws[1]:
        c.font = HEAD
        c.alignment = Alignment(vertical="top", wrap_text=True)
    for row in ws.iter_rows(min_row=2):
        for c in row:
            c.font = BODY
    ws.freeze_panes = "A2"
    ws.auto_filter.ref = ws.dimensions
    for i in range(1, ws.max_column + 1):
        w = max(len(str(ws.cell(r, i).value or ""))
                for r in range(1, min(ws.max_row, 200) + 1))
        ws.column_dimensions[get_column_letter(i)].width = min(max(w + 2, 10), 46)
    ws.row_dimensions[1].height = 30


def main():
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    wb = Workbook()
    toc = wb.active
    toc.title = "Contents"
    toc.append(["Supplementary Material 2. Supplementary Tables S1 to S9"])
    toc["A1"].font = TITLE
    toc.append([])
    toc.append(["Nam SW. How much of a tumor expression signature survives "
                "adjustment for the dominant tumor-adjacent axis?"])
    toc.append(["Every worksheet below is a verbatim copy of the named file in "
                "the deposited result archive (Zenodo, DOI 10.5281/zenodo.22658669)."])
    toc.append(["No value in this workbook was computed here."])
    toc.append([])
    toc.append(["Worksheet", "Data rows", "Columns", "Source file in the archive",
                "Contents"])

    made = []
    for name, fn, desc in SHEETS:
        with open(os.path.join(SRC, fn), newline="", encoding="utf-8") as f:
            rows = list(csv.reader(f))
        ws = wb.create_sheet(name)
        write_sheet(ws, rows)
        toc.append([name, len(rows) - 1, len(rows[0]), fn, desc])
        made.append((name, len(rows) - 1))

    for c in toc[7]:
        c.font = HEAD
    for row in toc.iter_rows(min_row=8):
        for c in row:
            c.font = BODY
            c.alignment = Alignment(vertical="top", wrap_text=True)
    for row in toc.iter_rows(min_row=1, max_row=6):
        for c in row:
            if c.font.sz != 12:
                c.font = BODY
    for col, w in zip("ABCDE", (22, 11, 9, 52, 76)):
        toc.column_dimensions[col].width = w
    toc.freeze_panes = "A8"

    wb.save(OUT)
    print(f"{os.path.basename(OUT)}: {len(made)} worksheets plus contents, "
          f"{sum(n for _, n in made)} data rows, "
          f"{os.path.getsize(OUT)/1e6:.2f} MB")


if __name__ == "__main__":
    main()
