"""io_utils.py -- loaders for GEO series matrices and count tables."""

from __future__ import annotations

import gzip
import io
import os
import re
import tarfile

import numpy as np
import pandas as pd


def _open_maybe_gzip(path: str):
    return gzip.open(path, "rt", errors="replace") if path.endswith(".gz") \
        else open(path, "rt", errors="replace")


def read_series_matrix(path: str) -> pd.DataFrame:
    """Parse the !Sample_* header block of a GEO series matrix file.

    Returns a DataFrame indexed by GSM accession. Each `characteristics_ch1`
    entry of the form "key: value" becomes its own column.
    """
    rows: dict[str, list[str]] = {}
    char_blocks: list[list[str]] = []
    with _open_maybe_gzip(path) as fh:
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                break
            if not line.startswith("!Sample_"):
                continue
            key, *vals = line.rstrip("\n").split("\t")
            key = key[len("!Sample_"):]
            vals = [v.strip().strip('"') for v in vals]
            if key == "characteristics_ch1":
                char_blocks.append(vals)
            else:
                rows.setdefault(key, vals)

    if "geo_accession" not in rows:
        raise ValueError(f"No !Sample_geo_accession line found in {path}")

    df = pd.DataFrame({k: v for k, v in rows.items()
                       if len(v) == len(rows["geo_accession"])})
    df = df.set_index("geo_accession")

    # explode characteristics into named columns
    for block in char_blocks:
        parsed = [(b.split(":", 1) + [""])[:2] if ":" in b else ("", b)
                  for b in block]
        keys = {k.strip().lower() for k, _ in parsed if k.strip()}
        if len(keys) == 1:
            col = keys.pop()
            df[col] = [v.strip() for _, v in parsed]
        else:  # heterogeneous block -- keep each cell as raw text
            df[f"characteristics_{len(df.columns)}"] = block
    df.columns = [re.sub(r"\s+", "_", c.strip().lower()) for c in df.columns]
    return df


def read_counts(path: str, sep: str | None = None) -> pd.DataFrame:
    """Read a counts matrix (genes x samples) from a flat file, .gz, or a
    RAW .tar of per-sample two-column files."""
    if path.endswith(".tar") or path.endswith(".tar.gz"):
        return _read_counts_from_tar(path)
    if sep is None:
        sep = "," if path.replace(".gz", "").endswith(".csv") else "\t"
    df = pd.read_csv(path, sep=sep, index_col=0, low_memory=False)
    df = df.loc[:, [c for c in df.columns
                    if pd.api.types.is_numeric_dtype(df[c])]]
    return df


def _read_counts_from_tar(path: str) -> pd.DataFrame:
    series: dict[str, pd.Series] = {}
    with tarfile.open(path) as tf:
        for member in tf.getmembers():
            if not member.isfile():
                continue
            raw = tf.extractfile(member).read()
            if member.name.endswith(".gz"):
                raw = gzip.decompress(raw)
            text = raw.decode("utf-8", errors="replace")
            sub = pd.read_csv(io.StringIO(text), sep="\t", header=None,
                              comment="#", engine="python")
            if sub.shape[1] < 2:
                continue
            # drop a header row if the second column is not numeric
            if not str(sub.iloc[0, 1]).replace(".", "", 1).isdigit():
                sub = sub.iloc[1:]
            s = pd.Series(pd.to_numeric(sub.iloc[:, 1], errors="coerce").values,
                          index=sub.iloc[:, 0].astype(str).values)
            gsm = re.match(r"(GSM\d+)", os.path.basename(member.name))
            series[gsm.group(1) if gsm else os.path.basename(member.name)] = s
    if not series:
        raise ValueError(f"No per-sample count files found inside {path}")
    return pd.DataFrame(series)


def _read_annotation(path: str) -> pd.DataFrame:
    """Read an identifier map, with or without a header.

    Accepts either a bare two-column file (id, symbol) or a table with
    named columns, such as NCBI's Human.GRCh38.p13.annot.tsv.gz, which
    carries GeneID, Symbol and EnsemblGeneID.
    """
    head = pd.read_csv(path, sep="\t", nrows=1, dtype=str)
    cols = {c.lower(): c for c in head.columns}
    has_header = any(k in cols for k in ("symbol", "geneid", "gene_name"))
    if not has_header:
        return pd.read_csv(path, sep="\t", header=None, dtype=str,
                           names=["gene_id", "symbol"])
    df = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
    low = {c.lower(): c for c in df.columns}
    sym = low.get("symbol") or low.get("gene_name")
    if sym is None:
        raise ValueError(f"No symbol column in {path}; saw {list(df.columns)}")
    frames = []
    for key in ("geneid", "ensemblgeneid", "gene_id", "ensembl_gene_id"):
        if key in low:
            sub = df[[low[key], sym]].copy()
            sub.columns = ["gene_id", "symbol"]
            # a cell may list several ids separated by commas
            sub = sub.assign(gene_id=sub["gene_id"].str.split(",")).explode("gene_id")
            frames.append(sub)
    if not frames:
        raise ValueError(f"No id column in {path}; saw {list(df.columns)}")
    out = pd.concat(frames, ignore_index=True).dropna()
    out["gene_id"] = out["gene_id"].str.strip()
    return out[out["gene_id"] != ""]


def to_symbols(counts: pd.DataFrame,
               annotation: str | None = None) -> pd.DataFrame:
    """Ensure the index is HGNC symbols.

    Handles three index flavors: symbols already, Ensembl gene IDs
    (ENSG..., with or without a version suffix), and Entrez GeneIDs
    (bare integers). The latter two need `annotation`.

    A suitable map for both is NCBI's Human.GRCh38.p13.annot.tsv.gz,
    downloadable from the GEO RNA-seq counts endpoint; or build one from
    the study's GTF:

        zcat gencode.v39.annotation.gtf.gz | awk '$3=="gene"' \
          | sed -n 's/.*gene_id "\([^"]*\)".*gene_name "\([^"]*\)".*/\1\t\2/p' \
          > ensembl_to_symbol.tsv
    """
    idx = counts.index.astype(str)
    is_ensembl = idx.str.match(r"^ENSG\d+").mean() > 0.5
    is_entrez = idx.str.fullmatch(r"\d+").mean() > 0.5

    if not (is_ensembl or is_entrez):
        out = counts.copy()
        if annotation is not None:
            amap = _read_annotation(annotation)
            keys = set(amap["gene_id"])
            covered = float(np.mean([g in keys for g in idx]))
            if covered > 0.001:
                # a platform-specific id map, e.g. Affymetrix probeset -> symbol
                m = dict(zip(amap["gene_id"], amap["symbol"]))
                mapped = [m.get(g) for g in idx]
                kept = sum(v is not None for v in mapped)
                print(f"  mapped {kept}/{len(mapped)} platform identifiers to "
                      f"symbols ({len(mapped) - kept} unmapped, dropped)")
                out.index = pd.Index([v if v is not None else "" for v in mapped])
                out = out.loc[out.index != ""]
                out.index = out.index.str.upper()
                return out.groupby(level=0).max()
        out.index = idx.str.upper()
        out = out.groupby(level=0).max()
        if annotation is not None:
            out = resolve_symbol_aliases(out, annotation)
        return out

    kind = "Ensembl gene IDs" if is_ensembl else "Entrez GeneIDs"
    if annotation is None:
        raise ValueError(
            f"Index looks like {kind} but no --annotation was given. "
            "Provide an id -> symbol map (see docstring).")

    amap = _read_annotation(annotation)
    keys = amap["gene_id"].str.split(".").str[0] if is_ensembl \
        else amap["gene_id"]
    m = dict(zip(keys, amap["symbol"]))
    stripped = idx.str.split(".").str[0] if is_ensembl else idx
    mapped = [m.get(g) for g in stripped]
    n_unmapped = sum(v is None for v in mapped)
    print(f"  mapped {len(mapped) - n_unmapped}/{len(mapped)} {kind} "
          f"to symbols ({n_unmapped} unmapped, dropped)")
    out = counts.copy()
    out.index = pd.Index([v if v is not None else "" for v in mapped])
    out = out.loc[out.index != ""]
    out.index = out.index.str.upper()
    return out.groupby(level=0).max()


SUMMARY_PREFIXES = ("__", "N_")   # HTSeq "__no_feature" etc.; STAR "N_unmapped"


def drop_summary_rows(counts: pd.DataFrame) -> pd.DataFrame:
    """Remove aligner summary rows.

    GEO per-sample count files from HTSeq end with __no_feature,
    __ambiguous, __too_low_aQual, __not_aligned and
    __alignment_not_unique; STAR's ReadsPerGene files start with N_*.
    These are read tallies, not genes. Leaving them in inflates every
    library size and therefore corrupts every CPM.
    """
    idx = counts.index.astype(str)
    mask = ~idx.str.startswith(SUMMARY_PREFIXES)
    dropped = int((~mask).sum())
    if dropped:
        print(f"  dropped {dropped} aligner summary rows "
              f"({', '.join(idx[~mask][:5])}...)")
    return counts.loc[mask]


def resolve_symbol_aliases(mat: pd.DataFrame,
                           annotation: str) -> pd.DataFrame:
    """Rename legacy gene symbols to their current HGNC symbol.

    Older GEO submissions use retired symbols -- MARC1/MARC2 for
    MTARC1/MTARC2, for example. A locked gene list written in current symbols
    silently loses those genes, which quietly degrades a module instead of
    failing. This resolves them from the Synonyms column of NCBI's annotation
    table.

    Only unambiguous aliases are applied, and only when the current symbol is
    not already present in the matrix, so nothing is ever overwritten.
    """
    df = pd.read_csv(annotation, sep="\t", dtype=str, low_memory=False)
    low = {c.lower(): c for c in df.columns}
    if "symbol" not in low or "synonyms" not in low:
        return mat
    cur = df[low["symbol"]].str.upper()
    syn = df[low["synonyms"]].fillna("")
    alias: dict[str, set] = {}
    for symbol, names in zip(cur, syn):
        for a in str(names).split("|"):
            a = a.strip().upper()
            if a and a != symbol:
                alias.setdefault(a, set()).add(symbol)

    present = set(mat.index)
    valid = set(cur)
    renames = {}
    for g in mat.index:
        if g in valid:
            continue                       # already a current symbol
        targets = alias.get(g)
        if targets and len(targets) == 1:
            t = next(iter(targets))
            if t not in present:           # never overwrite a real row
                renames[g] = t
    if not renames:
        return mat
    shown = ", ".join(f"{k}->{v}" for k, v in list(renames.items())[:5])
    print(f"  resolved {len(renames)} legacy gene symbols ({shown}"
          f"{', ...' if len(renames) > 5 else ''})")
    out = mat.rename(index=renames)
    return out.groupby(level=0).max()


def log2_cpm(counts: pd.DataFrame, min_cpm: float = 1.0,
             min_samples_frac: float = 0.20) -> pd.DataFrame:
    """Filter lowly expressed genes, then return log2(CPM + 1)."""
    counts = drop_summary_rows(counts).fillna(0)
    lib = counts.sum(axis=0)
    if (lib <= 0).any():
        bad = list(lib.index[lib <= 0])
        raise ValueError(f"Samples with zero library size: {bad}")
    cpm = counts.divide(lib, axis=1) * 1e6
    keep = (cpm >= min_cpm).sum(axis=1) >= max(
        3, int(round(min_samples_frac * counts.shape[1])))
    return np.log2(cpm.loc[keep] + 1.0)


def quantile_normalize(mat: pd.DataFrame) -> pd.DataFrame:
    """Force every sample to share one distribution.

    GSE164760's probeset table is on a linear scale and is NOT equalized
    across samples (total-intensity CV about 8%). Quantile normalization is
    the standard remedy for array data in that state, and it is applied
    before any index is computed so that between-group contrasts are not
    driven by per-array scaling.
    """
    ranks = mat.rank(axis=0, method="average")
    reference = np.sort(mat.to_numpy(), axis=0).mean(axis=1)
    idx = np.clip(np.round(ranks.to_numpy()).astype(int) - 1,
                  0, len(reference) - 1)
    return pd.DataFrame(reference[idx], index=mat.index, columns=mat.columns)


def already_logged(mat: pd.DataFrame) -> bool:
    """Heuristic: microarray/RMA matrices arrive already log-scaled."""
    return float(np.nanmax(mat.values)) < 40.0
