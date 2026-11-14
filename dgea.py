import csv
import gzip
import re
from io import StringIO

import numpy as np
import pandas as pd
from scipy.stats import ttest_rel


def benjamini_hochberg(p_values):
    """Adjust a collection of p-values while preserving missing values."""
    p = np.asarray(p_values, dtype=float)
    result = np.full(p.shape, np.nan)
    valid = np.isfinite(p)
    vals = p[valid]
    if not len(vals):
        return result

    order = np.argsort(vals)
    ranked = vals[order] * len(vals) / np.arange(1, len(vals) + 1)
    adjusted = np.minimum.accumulate(ranked[::-1])[::-1]
    back = np.empty_like(adjusted)
    back[order] = np.minimum(adjusted, 1.0)
    result[valid] = back
    return result


def _csv_fields(line):
    return next(csv.reader([line], delimiter="\t", quotechar='"'))[1:]


def read_series_matrix(path):
    """Read a GEO series matrix and the sample fields used by this project."""
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as handle:
        lines = handle.readlines()

    titles = next(line for line in lines if line.startswith("!Sample_title\t"))
    accessions = next(line for line in lines if line.startswith("!Sample_geo_accession\t"))
    sources = next(line for line in lines if line.startswith("!Sample_source_name_ch1\t"))
    processing = next(
        (line for line in lines if line.startswith("!Sample_data_processing\t")),
        None,
    )

    title_values = _csv_fields(titles)
    gsm_values = _csv_fields(accessions)
    source_values = _csv_fields(sources)
    processing_values = _csv_fields(processing) if processing else [""] * len(gsm_values)

    start = next(i for i, line in enumerate(lines) if line.startswith('"ID_REF"\t'))
    end = next(
        (i for i in range(start + 1, len(lines)) if lines[i].startswith("!series_matrix_table_end")),
        len(lines),
    )
    matrix = pd.read_csv(
        StringIO("".join(lines[start:end])),
        sep="\t",
        index_col=0,
        low_memory=False,
        quotechar='"',
    )
    metadata = pd.DataFrame(
        {
            "gsm": gsm_values,
            "title": title_values,
            "source": source_values,
            "processing": processing_values,
        }
    )
    return matrix.apply(pd.to_numeric, errors="coerce"), metadata


def read_platform_annotation(path):
    """Read GEO's GPL96 annotation and keep unambiguous gene mappings."""
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as handle:
        lines = handle.readlines()

    start = next(i for i, line in enumerate(lines) if line.startswith("ID\tGene title\tGene symbol\t"))
    end = next(
        (i for i in range(start + 1, len(lines)) if lines[i].startswith("!platform_table_end")),
        len(lines),
    )
    annotation = pd.read_csv(
        StringIO("".join(lines[start:end])),
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )[["ID", "Gene title", "Gene symbol"]]
    annotation.columns = ["probe_id", "gene_title", "gene_symbol"]

    symbol = annotation.gene_symbol.str.strip()
    title = annotation.gene_title.str.strip()
    unambiguous = symbol.ne("") & symbol.ne("---") & ~symbol.str.contains("///", regex=False)
    annotation.loc[~unambiguous, "gene_symbol"] = pd.NA
    annotation.loc[title.isin(["", "---"]), "gene_title"] = pd.NA
    return annotation.drop_duplicates("probe_id")


def collapse_probes_to_genes(matrix, annotation):
    """Average probes that map unambiguously to the same gene symbol."""
    mapped = annotation.set_index("probe_id").reindex(matrix.index)
    keep = mapped.gene_symbol.notna()
    if not keep.any():
        raise ValueError("No matrix probes mapped to unambiguous gene symbols")

    expression = matrix.loc[keep].copy()
    expression.insert(0, "gene_symbol", mapped.loc[keep, "gene_symbol"].to_numpy())
    gene_matrix = expression.groupby("gene_symbol", sort=True)[matrix.columns].mean()

    titles = (
        mapped.loc[keep, ["gene_symbol", "gene_title"]]
        .dropna()
        .drop_duplicates("gene_symbol")
        .set_index("gene_symbol")
        .gene_title
        .reindex(gene_matrix.index)
    )
    return gene_matrix, titles


def log2_transform(matrix):
    """Log2-transform positive normalized intensities."""
    values = matrix.to_numpy(dtype=float)
    finite = values[np.isfinite(values)]
    if finite.size == 0 or (finite <= 0).any():
        raise ValueError("Log2 transformation requires positive expression values")
    return np.log2(matrix.astype(float))


def classify_gse10072(metadata):
    """Keep subject IDs with exactly one tumor and one normal profile."""
    rows = []
    for row in metadata.itertuples(index=False):
        source = row.source.lower()
        if "adenocarcinoma" in source:
            group = "tumor"
        elif "normal lung" in source:
            group = "normal"
        else:
            continue
        subject = re.search(r"GT\d+", row.title)
        if subject:
            rows.append((row.gsm, subject.group(0), group))

    frame = pd.DataFrame(rows, columns=["gsm", "pair", "group"])
    counts = frame.groupby(["pair", "group"]).size().unstack(fill_value=0)
    keep = counts[
        (counts.get("normal", 0) == 1) & (counts.get("tumor", 0) == 1)
    ].index
    return frame[frame.pair.isin(keep)].copy()


def classify_gse7670(metadata):
    """Keep only unique numbered tumor/adjacent-normal patient pairs."""
    rows = []
    for row in metadata.itertuples(index=False):
        match = re.match(r"^(\d+)(N|T):", row.title)
        if match:
            rows.append(
                (row.gsm, match.group(1), "normal" if match.group(2) == "N" else "tumor")
            )
    frame = pd.DataFrame(rows, columns=["gsm", "pair", "group"])
    counts = frame.groupby(["pair", "group"]).size().unstack(fill_value=0)
    keep = counts[(counts.get("normal", 0) == 1) & (counts.get("tumor", 0) == 1)].index
    return frame[frame.pair.isin(keep)].copy()


def _results_table(matrix, normal, tumor, statistic, p_values, gene_titles):
    difference = tumor.mean(axis=1) - normal.mean(axis=1)
    difference_sd = (
        (tumor.to_numpy() - normal.to_numpy()).std(axis=1, ddof=1)
        if normal.shape == tumor.shape
        else None
    )
    result = pd.DataFrame(
        {
            "gene_symbol": matrix.index,
            "gene_title": gene_titles.reindex(matrix.index).to_numpy(),
            "mean_log2_expression_normal": normal.mean(axis=1).to_numpy(),
            "mean_log2_expression_tumor": tumor.mean(axis=1).to_numpy(),
            "log2_fold_change_tumor_vs_normal": difference.to_numpy(),
            "t_statistic": statistic,
            "p_value": p_values,
            "p_adj_bh": benjamini_hochberg(p_values),
        }
    )
    if difference_sd is not None:
        result.insert(
            5,
            "cohens_dz",
            np.divide(
                difference.to_numpy(),
                difference_sd,
                out=np.full(len(difference), np.nan),
                where=difference_sd != 0,
            ),
        )
    return result.sort_values(["p_adj_bh", "p_value"], na_position="last").reset_index(drop=True)


def analyze_paired(matrix, samples, gene_titles):
    """Run a two-sided paired t-test for matched tumor/normal samples."""
    wide = samples.pivot(index="pair", columns="group", values="gsm").sort_index()
    if wide[["normal", "tumor"]].isna().any().any():
        raise ValueError("Paired analysis requires one normal and one tumor sample per pair")
    normal = matrix[wide.normal]
    tumor = matrix[wide.tumor]
    statistic, p_values = ttest_rel(tumor, normal, axis=1, nan_policy="omit")
    return _results_table(matrix, normal, tumor, statistic, p_values, gene_titles)
