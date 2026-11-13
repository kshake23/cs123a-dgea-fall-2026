import gzip
import re

import numpy as np
import pandas as pd
from scipy.stats import ttest_rel


def benjamini_hochberg(p_values):
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


def read_series_matrix(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as handle:
        lines = handle.readlines()
    titles = next(line for line in lines if line.startswith("!Sample_title\t"))
    accessions = next(line for line in lines if line.startswith("!Sample_geo_accession\t"))
    sources = next(line for line in lines if line.startswith("!Sample_source_name_ch1\t"))
    def fields(line):
        return next(csv_reader([line]))[1:]
    title_values, gsm_values = fields(titles), fields(accessions)
    start = next(i for i, line in enumerate(lines) if line.startswith('"ID_REF"\t'))
    end = next((i for i in range(start + 1, len(lines)) if lines[i].startswith("!series_matrix_table_end")), len(lines))
    rows = lines[start + 1:end]
    from io import StringIO
    matrix = pd.read_csv(StringIO("".join([lines[start], *rows])), sep="\t", index_col=0, low_memory=False, quotechar='"')
    matrix.columns = gsm_values
    source_values = fields(sources)
    return matrix.apply(pd.to_numeric, errors="coerce"), pd.DataFrame({"gsm": gsm_values, "title": title_values, "source": source_values})


def csv_reader(lines):
    import csv
    return csv.reader(lines, delimiter="\t", quotechar='"')


def classify_gse10072(metadata):
    rows = []
    for row in metadata.itertuples(index=False):
        source = row.source
        value = source.lower()
        group = "tumor" if "adenocarcinoma" in value else "normal" if "normal lung" in value else None
        match = re.search(r"GT\d+", row.title)
        if group and match:
            rows.append((row.gsm, match.group(0), group))
    frame = pd.DataFrame(rows, columns=["gsm", "pair", "group"])
    counts = frame.groupby(["pair", "group"]).size().unstack(fill_value=0)
    complete = counts[(counts.get("normal", 0) == 1) & (counts.get("tumor", 0) == 1)].index
    return frame[frame.pair.isin(complete)]


def classify_gse7670(metadata):
    rows = []
    for row in metadata.itertuples(index=False):
        match = re.match(r"^(\d+)(N|T):", row.title)
        if match:
            rows.append((row.gsm, match.group(1), "normal" if match.group(2) == "N" else "tumor"))
    frame = pd.DataFrame(rows, columns=["gsm", "pair", "group"])
    counts = frame.groupby(["pair", "group"]).size().unstack(fill_value=0)
    keep = counts[(counts.get("normal", 0) == 1) & (counts.get("tumor", 0) == 1)].index
    return frame[frame.pair.isin(keep)]


def analyze_paired(matrix, samples):
    wide = samples.pivot(index="pair", columns="group", values="gsm")
    normal, tumor = matrix[wide.normal], matrix[wide.tumor]
    stat, p = ttest_rel(tumor, normal, axis=1, nan_policy="omit")
    return results_table(matrix, tumor.mean(axis=1) - normal.mean(axis=1), stat, p)


def results_table(matrix, difference, statistic, p_values):
    return pd.DataFrame({"probe_id": matrix.index, "mean_log_expression_difference_tumor_minus_normal": difference,
                         "t_statistic": statistic, "p_value": p_values,
                         "p_adj_bh": benjamini_hochberg(p_values)}).sort_values("p_adj_bh")
