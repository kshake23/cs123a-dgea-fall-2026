from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from dgea import (
    classify_gse7670,
    collapse_probes_to_genes,
    log2_transform,
    read_platform_annotation,
    read_series_matrix,
)


DATA = ROOT / "data"
RESULTS = ROOT / "results"
FC_CUTOFF = 1.0
FDR_CUTOFF = 0.05
COLORS = {
    "background": "#f7f5f0",
    "text": "#17253d",
    "muted": "#9aa3ad",
    "blue": "#2878b5",
    "orange": "#e07a3f",
}


def _style_axes(ax):
    ax.set_facecolor(COLORS["background"])
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(colors=COLORS["text"])
    ax.xaxis.label.set_color(COLORS["text"])
    ax.yaxis.label.set_color(COLORS["text"])
    ax.title.set_color(COLORS["text"])


def make_volcano(accession):
    table = pd.read_csv(RESULTS / f"{accession}_dgea.csv")
    x = table["log2_fold_change_tumor_vs_normal"]
    y = -np.log10(table["p_adj_bh"].clip(lower=np.finfo(float).tiny))
    highlighted = (table["p_adj_bh"] < FDR_CUTOFF) & (x.abs() >= FC_CUTOFF)

    fig, ax = plt.subplots(figsize=(8, 4.6), facecolor=COLORS["background"])
    _style_axes(ax)
    ax.scatter(x[~highlighted], y[~highlighted], s=7, alpha=0.25, color=COLORS["muted"])
    ax.scatter(
        x[highlighted],
        y[highlighted],
        s=8,
        alpha=0.55,
        color=COLORS["orange"],
        label="FDR < .05; |log2FC| >= 1",
    )
    ax.axhline(-np.log10(FDR_CUTOFF), color=COLORS["text"], linestyle="--", linewidth=1)
    ax.axvline(-FC_CUTOFF, color=COLORS["text"], linestyle=":", linewidth=1)
    ax.axvline(FC_CUTOFF, color=COLORS["text"], linestyle=":", linewidth=1)
    ax.set_xlabel("log2 fold change (tumor - normal)")
    ax.set_ylabel("-log10(BH-adjusted p-value)")
    ax.set_title(f"{accession} differential expression", loc="left", weight="bold")
    ax.legend(frameon=False, markerscale=2, fontsize=9, loc="upper right", borderaxespad=0.4)
    fig.tight_layout()
    fig.savefig(RESULTS / f"{accession}_volcano.png", dpi=180, facecolor=fig.get_facecolor())
    plt.close(fig)


def make_concordance():
    first = pd.read_csv(RESULTS / "GSE10072_dgea.csv").set_index("gene_symbol")
    second = pd.read_csv(RESULTS / "GSE7670_dgea.csv").set_index("gene_symbol")
    first = first[["log2_fold_change_tumor_vs_normal", "p_adj_bh"]].add_suffix("_10072")
    second = second[["log2_fold_change_tumor_vs_normal", "p_adj_bh"]].add_suffix("_7670")
    joined = first.join(second, how="inner")
    x = joined["log2_fold_change_tumor_vs_normal_10072"]
    y = joined["log2_fold_change_tumor_vs_normal_7670"]
    both = (joined["p_adj_bh_10072"] < FDR_CUTOFF) & (
        joined["p_adj_bh_7670"] < FDR_CUTOFF
    )
    same_direction = np.sign(x) == np.sign(y)
    concordant = both & same_direction
    rho, _ = spearmanr(x, y)
    same_direction_count = int(concordant.sum())
    print(
        f"Across {len(joined):,} shared genes: Spearman rho = {rho:.3f}; "
        f"{both.sum():,} significant in both studies; "
        f"{same_direction_count:,} significant in both with the same direction."
    )

    fig, ax = plt.subplots(figsize=(7, 5.4), facecolor=COLORS["background"])
    _style_axes(ax)
    ax.scatter(x[~concordant], y[~concordant], s=7, alpha=0.18, color=COLORS["muted"])
    ax.scatter(
        x[concordant],
        y[concordant],
        s=8,
        alpha=0.38,
        color=COLORS["blue"],
        label="FDR < 0.05 in both, same direction",
    )
    low = min(x.min(), y.min())
    high = max(x.max(), y.max())
    ax.plot([low, high], [low, high], color=COLORS["text"], linewidth=1, linestyle="--")
    ax.axhline(0, color=COLORS["text"], linewidth=0.7, alpha=0.5)
    ax.axvline(0, color=COLORS["text"], linewidth=0.7, alpha=0.5)

    label_offsets = {"FABP4": (-10, 9), "SFTPC": (12, 0)}
    for gene in ("SPP1", "MMP1", "COL11A1", "AGER", "FABP4", "SFTPC"):
        if gene in joined.index:
            offset = label_offsets.get(gene, (4, 4))
            ax.annotate(
                gene,
                (x.loc[gene], y.loc[gene]),
                xytext=offset,
                textcoords="offset points",
                ha="right" if offset[0] < 0 else "left",
                arrowprops=(
                    {"arrowstyle": "-", "color": COLORS["text"], "lw": 0.6}
                    if gene in label_offsets
                    else None
                ),
                fontsize=8,
                color=COLORS["text"],
            )

    ax.text(
        0.03,
        0.96,
        f"Spearman rho = {rho:.3f}\n{same_direction_count:,} genes agree in direction",
        transform=ax.transAxes,
        va="top",
        color=COLORS["text"],
        fontsize=10,
    )
    ax.set_xlabel("GSE10072 log2 fold change")
    ax.set_ylabel("GSE7670 log2 fold change")
    ax.set_title("Fold-change agreement across datasets", loc="left", weight="bold")
    ax.legend(frameon=False, markerscale=2, loc="lower right")
    fig.tight_layout()
    fig.savefig(RESULTS / "cross_dataset_concordance.png", dpi=180, facecolor=fig.get_facecolor())
    plt.close(fig)


def make_gse7670_paired_plot():
    matrix, metadata = read_series_matrix(DATA / "GSE7670_series_matrix.txt.gz")
    annotation = read_platform_annotation(DATA / "GPL96.annot.gz")
    samples = classify_gse7670(metadata)
    gene_matrix, _ = collapse_probes_to_genes(log2_transform(matrix[samples.gsm]), annotation)
    wide = samples.pivot(index="pair", columns="group", values="gsm").sort_index()

    fig, axes = plt.subplots(1, 2, figsize=(7, 4.6), facecolor=COLORS["background"], sharey=False)
    for ax, gene, color in zip(axes, ("AGER", "SPP1"), (COLORS["blue"], COLORS["orange"])):
        normal = gene_matrix.loc[gene, wide["normal"]].to_numpy(dtype=float)
        tumor = gene_matrix.loc[gene, wide["tumor"]].to_numpy(dtype=float)
        _style_axes(ax)
        for before, after in zip(normal, tumor):
            ax.plot([0, 1], [before, after], color=color, alpha=0.24, linewidth=1)
        ax.scatter(np.zeros_like(normal), normal, color=color, s=18, alpha=0.65)
        ax.scatter(np.ones_like(tumor), tumor, color=color, s=18, alpha=0.65)
        ax.plot(
            [0, 1],
            [normal.mean(), tumor.mean()],
            color=COLORS["text"],
            linewidth=3,
            marker="o",
            markersize=6,
        )
        ax.set_xticks([0, 1], ["Normal", "Tumor"])
        ax.set_title(gene, loc="left", weight="bold")
        ax.set_ylabel("log2 expression")

    fig.suptitle("GSE7670 matched-patient expression", x=0.08, ha="left", color=COLORS["text"], weight="bold")
    fig.text(0.08, 0.02, "Each thin line is one of 27 matched patients; the thick line connects group means.", color=COLORS["text"], fontsize=9)
    fig.tight_layout(rect=(0, 0.06, 1, 0.94))
    fig.savefig(RESULTS / "GSE7670_paired_genes.png", dpi=180, facecolor=fig.get_facecolor())
    plt.close(fig)


if __name__ == "__main__":
    for accession in ("GSE10072", "GSE7670"):
        make_volcano(accession)
    make_concordance()
    make_gse7670_paired_plot()
