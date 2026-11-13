from pathlib import Path
import gzip
import urllib.request
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from dgea import analyze_paired, classify_gse10072, classify_gse7670, read_series_matrix
DATA = ROOT / "data"
RESULTS = ROOT / "results"
SOURCES = {
    "GSE10072": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE10nnn/GSE10072/matrix/GSE10072_series_matrix.txt.gz",
    "GSE7670": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE7nnn/GSE7670/matrix/GSE7670_series_matrix.txt.gz",
}


def download(accession, url):
    path = DATA / f"{accession}_series_matrix.txt.gz"
    if not path.exists():
        urllib.request.urlretrieve(url, path)
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        if f'!Series_geo_accession\t"{accession}"' not in "".join(handle.readline() for _ in range(3)):
            raise ValueError(f"{path} is not the expected GEO series matrix")
    return path


def plot_volcano(table, accession):
    x = table.mean_log_expression_difference_tumor_minus_normal
    y = -np.log10(table.p_adj_bh.clip(lower=np.finfo(float).tiny))
    plt.figure(figsize=(7, 5))
    plt.scatter(x, y, s=5, alpha=.45)
    plt.axhline(-np.log10(.05), color="gray", linestyle="--")
    plt.xlabel("Mean log-expression difference (tumor − normal)")
    plt.ylabel("−log10(BH adjusted p-value)")
    plt.title(f"{accession} differential expression")
    plt.tight_layout()
    plt.savefig(RESULTS / f"{accession}_volcano.png", dpi=160)
    plt.close()


def main():
    DATA.mkdir(exist_ok=True)
    RESULTS.mkdir(exist_ok=True)
    for accession, url in SOURCES.items():
        matrix, metadata = read_series_matrix(download(accession, url))
        if accession == "GSE10072":
            samples = classify_gse10072(metadata)
            matrix = matrix[samples.gsm]
            table = analyze_paired(matrix, samples)
            manifest = metadata.merge(samples, on="gsm", how="inner")
            print(f"GSE10072: retaining {manifest.pair.nunique()} complete pairs; excluding {len(metadata) - len(manifest)} unmatched/unclassified samples")
        else:
            samples = classify_gse7670(metadata)
            matrix = matrix[samples.gsm]
            table = analyze_paired(matrix, samples)
            manifest = metadata.merge(samples, on="gsm", how="inner")
            if manifest.pair.nunique() != 27:
                raise ValueError(f"Expected 27 matched pairs; found {manifest.pair.nunique()}")
        table.to_csv(RESULTS / f"{accession}_dgea.csv", index=False)
        manifest.to_csv(RESULTS / f"{accession}_samples.csv", index=False)
        plot_volcano(table, accession)
        print(f"{accession}: {len(table)} probes; {len(manifest)} samples; significant FDR<0.05: {(table.p_adj_bh < .05).sum()}")


if __name__ == "__main__":
    main()
