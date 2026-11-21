from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"


def make_volcano(accession):
    table = pd.read_csv(RESULTS / f"{accession}_dgea.csv")
    x = table["log2_fold_change_tumor_vs_normal"]
    y = -np.log10(table["p_adj_bh"].clip(lower=np.finfo(float).tiny))
    significant = table["p_adj_bh"] < 0.05

    plt.figure(figsize=(7, 5))
    plt.scatter(x[~significant], y[~significant], s=5, alpha=0.25, label="Not significant")
    plt.scatter(x[significant], y[significant], s=5, alpha=0.35, label="BH FDR < 0.05")
    plt.axhline(-np.log10(0.05), color="gray", linestyle="--", linewidth=1)
    plt.xlabel("log2 fold change (tumor − normal)")
    plt.ylabel("−log10(BH-adjusted p-value)")
    plt.title(f"{accession} gene-level differential expression")
    plt.legend(frameon=False, markerscale=2)
    plt.tight_layout()
    plt.savefig(RESULTS / f"{accession}_volcano.png", dpi=180)
    plt.close()


if __name__ == "__main__":
    for accession in ("GSE10072", "GSE7670"):
        make_volcano(accession)
