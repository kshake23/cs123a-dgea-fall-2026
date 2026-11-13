from pathlib import Path
import gzip
import sys
import urllib.request

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from adjustText import adjust_text
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from dgea import (
    analyze_paired,
    analyze_unpaired,
    classify_gse10072,
    classify_gse7670,
    collapse_probes_to_genes,
    log2_transform,
    read_platform_annotation,
    read_series_matrix,
)

DATA = ROOT / "data"
RESULTS = ROOT / "results"
SOURCES = {
    "GSE10072": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE10nnn/GSE10072/matrix/GSE10072_series_matrix.txt.gz",
    "GSE7670": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE7nnn/GSE7670/matrix/GSE7670_series_matrix.txt.gz",
}
ANNOTATION_URL = "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL96/annot/GPL96.annot.gz"


def download_gzip(path, url, required_text):
    if not path.exists():
        urllib.request.urlretrieve(url, path)
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        header = "".join(handle.readline() for _ in range(35))
    if required_text not in header:
        raise ValueError(f"{path} failed its source validation check")
    return path


def style_axes(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="#E5E7EB", linewidth=0.7, alpha=0.8)
    ax.set_axisbelow(True)


def plot_volcano(table, accession):
    x = table.log2_fold_change_tumor_vs_normal
    y = -np.log10(table.p_adj_bh.clip(lower=np.finfo(float).tiny))
    significant = table.p_adj_bh.lt(0.05) & x.abs().ge(1)
    colors = np.where(
        significant & x.gt(0),
        "#C7472C",
        np.where(significant, "#2878B5", "#B8BDC5"),
    )

    fig, ax = plt.subplots(figsize=(8.4, 5.8))
    ax.scatter(x, y, s=14, c=colors, alpha=0.72, linewidths=0)
    ax.axhline(-np.log10(0.05), color="#565B63", linestyle="--", linewidth=1)
    ax.axvline(-1, color="#8B9098", linestyle=":", linewidth=1)
    ax.axvline(1, color="#8B9098", linestyle=":", linewidth=1)

    candidates = pd.concat(
        [
            table.loc[significant & x.lt(0)].nsmallest(5, "p_adj_bh"),
            table.loc[significant & x.gt(0)].nsmallest(4, "p_adj_bh"),
        ]
    )
    texts = []
    for row in candidates.itertuples(index=False):
        texts.append(
            ax.text(
                row.log2_fold_change_tumor_vs_normal,
                -np.log10(max(row.p_adj_bh, np.finfo(float).tiny)),
                row.gene_symbol,
                fontsize=8,
            )
        )
    adjust_text(
        texts,
        ax=ax,
        expand=(1.08, 1.18),
        arrowprops={"arrowstyle": "-", "color": "#7B8088", "lw": 0.55},
    )
    ax.set(
        xlabel="log2 fold change (tumor vs. normal)",
        ylabel="-log10(BH-adjusted p-value)",
        title=f"{accession}: differential gene expression",
    )
    style_axes(ax)
    fig.tight_layout()
    fig.savefig(RESULTS / f"{accession}_volcano.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


def plot_gse7670_paired_genes(matrix, samples, table):
    top = table.loc[
        table.p_adj_bh.lt(0.05)
        & table.log2_fold_change_tumor_vs_normal.abs().ge(1)
    ].head(6)
    wide = samples.pivot(index="pair", columns="group", values="gsm").sort_index()
    fig, axes = plt.subplots(2, 3, figsize=(12, 7.2), sharex=False)
    for ax, row in zip(axes.flat, top.itertuples(index=False)):
        normal = matrix.loc[row.gene_symbol, wide.normal].to_numpy()
        tumor = matrix.loc[row.gene_symbol, wide.tumor].to_numpy()
        for normal_value, tumor_value in zip(normal, tumor):
            ax.plot(
                [0, 1],
                [normal_value, tumor_value],
                color="#C9CDD3",
                linewidth=0.8,
                alpha=0.8,
            )
        ax.scatter(np.zeros(len(normal)), normal, s=18, color="#2878B5", zorder=3)
        ax.scatter(np.ones(len(tumor)), tumor, s=18, color="#C7472C", zorder=3)
        ax.set_xticks([0, 1], ["Adjacent\nnormal", "Tumor"])
        ax.set_title(
            f"{row.gene_symbol}  (log2FC {row.log2_fold_change_tumor_vs_normal:+.2f})",
            fontsize=10,
        )
        ax.set_ylabel("log2 expression")
        style_axes(ax)
    fig.suptitle(
        "GSE7670: paired-patient expression for leading genes",
        fontsize=15,
        y=1.01,
    )
    fig.tight_layout()
    fig.savefig(
        RESULTS / "GSE7670_paired_top_genes.png",
        dpi=220,
        bbox_inches="tight",
    )
    plt.close(fig)


def compare_datasets(tables):
    left = tables["GSE10072"][
        ["gene_symbol", "log2_fold_change_tumor_vs_normal", "p_adj_bh"]
    ].rename(
        columns={
            "log2_fold_change_tumor_vs_normal": "log2fc_GSE10072",
            "p_adj_bh": "fdr_GSE10072",
        }
    )
    right = tables["GSE7670"][
        ["gene_symbol", "gene_title", "log2_fold_change_tumor_vs_normal", "p_adj_bh"]
    ].rename(
        columns={
            "log2_fold_change_tumor_vs_normal": "log2fc_GSE7670",
            "p_adj_bh": "fdr_GSE7670",
        }
    )
    comparison = left.merge(right, on="gene_symbol", how="inner")
    comparison["significant_both"] = comparison.fdr_GSE10072.lt(
        0.05
    ) & comparison.fdr_GSE7670.lt(0.05)
    comparison["same_direction"] = np.sign(
        comparison.log2fc_GSE10072
    ) == np.sign(comparison.log2fc_GSE7670)
    comparison["mean_abs_log2fc"] = comparison[
        ["log2fc_GSE10072", "log2fc_GSE7670"]
    ].abs().mean(axis=1)
    comparison = comparison.sort_values(
        ["significant_both", "same_direction", "mean_abs_log2fc"],
        ascending=[False, False, False],
    )
    comparison.to_csv(RESULTS / "cross_dataset_comparison.csv", index=False)

    rho, rho_p = spearmanr(
        comparison.log2fc_GSE10072, comparison.log2fc_GSE7670
    )
    concordant = comparison.significant_both & comparison.same_direction
    fig, ax = plt.subplots(figsize=(7.2, 6.4))
    ax.scatter(
        comparison.log2fc_GSE10072,
        comparison.log2fc_GSE7670,
        s=11,
        color="#B8BDC5",
        alpha=0.45,
        linewidths=0,
        label="Shared genes",
    )
    ax.scatter(
        comparison.loc[concordant, "log2fc_GSE10072"],
        comparison.loc[concordant, "log2fc_GSE7670"],
        s=18,
        color="#6B4FA1",
        alpha=0.72,
        linewidths=0,
        label="FDR < 0.05 in both, same direction",
    )
    limit = float(
        np.nanmax(
            np.abs(
                comparison[["log2fc_GSE10072", "log2fc_GSE7670"]].to_numpy()
            )
        )
    )
    ax.plot(
        [-limit, limit],
        [-limit, limit],
        color="#565B63",
        linestyle="--",
        linewidth=1,
    )
    ax.axhline(0, color="#8B9098", linewidth=0.7)
    ax.axvline(0, color="#8B9098", linewidth=0.7)
    top = comparison.loc[concordant].nlargest(9, "mean_abs_log2fc")
    texts = []
    for row in top.itertuples(index=False):
        texts.append(
            ax.text(
                row.log2fc_GSE10072,
                row.log2fc_GSE7670,
                row.gene_symbol,
                fontsize=8,
            )
        )
    adjust_text(
        texts,
        ax=ax,
        expand=(1.08, 1.14),
        arrowprops={"arrowstyle": "-", "color": "#7B8088", "lw": 0.55},
    )
    ax.set(
        xlabel="GSE10072 log2 fold change",
        ylabel="GSE7670 log2 fold change",
        title=f"Cross-dataset effect agreement (Spearman rho = {rho:.2f})",
    )
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    style_axes(ax)
    fig.tight_layout()
    fig.savefig(
        RESULTS / "cross_dataset_concordance.png",
        dpi=220,
        bbox_inches="tight",
    )
    plt.close(fig)
    return comparison, rho, rho_p


def write_summary(tables, comparison, rho, rho_p):
    lines = [
        "# Analysis summary",
        "",
        "The two studies were analyzed separately and compared only after testing. GSE10072 used a two-sided Welch t-test on 58 tumor and 49 normal samples. GSE7670 used a two-sided paired t-test on 27 matched tumor/adjacent-normal patient pairs. BH adjustment controlled the false discovery rate.",
        "",
        "## Results",
        "",
    ]
    for accession, table in tables.items():
        fdr = table.p_adj_bh.lt(0.05)
        strong = fdr & table.log2_fold_change_tumor_vs_normal.abs().ge(1)
        up = strong & table.log2_fold_change_tumor_vs_normal.gt(0)
        down = strong & table.log2_fold_change_tumor_vs_normal.lt(0)
        lines.append(
            f"- **{accession}:** {len(table):,} unambiguously mapped genes tested; "
            f"{fdr.sum():,} at FDR < 0.05; {strong.sum():,} also had "
            f"|log2 fold change| >= 1 ({up.sum():,} higher and {down.sum():,} lower "
            "in tumor)."
        )
    both = comparison.significant_both
    concordant = both & comparison.same_direction
    lines.extend(
        [
            f"- **Cross-dataset:** {both.sum():,} genes were significant in both studies; "
            f"{concordant.sum():,} of those changed in the same direction.",
            f"- Across all {len(comparison):,} shared mapped genes, log2 fold changes had "
            f"Spearman rho = {rho:.3f} "
            f"(p {'< 1e-300' if rho_p == 0 else f'= {rho_p:.3g}'}).",
            "",
            "## Leading concordant genes",
            "",
            "| Gene | GSE10072 log2FC | GSE7670 log2FC | GSE10072 FDR | GSE7670 FDR |",
            "|---|---:|---:|---:|---:|",
        ]
    )
    top = comparison.loc[concordant].nlargest(12, "mean_abs_log2fc")
    for row in top.itertuples(index=False):
        lines.append(
            f"| {row.gene_symbol} | {row.log2fc_GSE10072:.2f} | "
            f"{row.log2fc_GSE7670:.2f} | {row.fdr_GSE10072:.2e} | "
            f"{row.fdr_GSE7670:.2e} |"
        )
    lines.extend(
        [
            "",
            "## Interpretation limits",
            "",
            "These are observational microarray results, not evidence that a gene causes lung cancer. The studies differ in cohort and preprocessing, and the GPL96 annotation is an older array definition. Only probes with a single current GEO gene-symbol mapping were retained; probes mapping to multiple genes were excluded, and multiple probes for the same gene were averaged before testing.",
            "",
        ]
    )
    (RESULTS / "analysis_summary.md").write_text(
        "\n".join(lines), encoding="utf-8"
    )


def main():
    DATA.mkdir(exist_ok=True)
    RESULTS.mkdir(exist_ok=True)
    annotation_path = download_gzip(
        DATA / "GPL96.annot.gz",
        ANNOTATION_URL,
        "!Annotation_platform = GPL96",
    )
    annotation = read_platform_annotation(annotation_path)

    tables = {}
    matrices = {}
    samples_by_accession = {}
    for accession, url in SOURCES.items():
        matrix_path = download_gzip(
            DATA / f"{accession}_series_matrix.txt.gz",
            url,
            f'!Series_geo_accession\t"{accession}"',
        )
        matrix, metadata = read_series_matrix(matrix_path)

        if accession == "GSE10072":
            samples = classify_gse10072(metadata)
            matrix = matrix[samples.gsm]
            if matrix.max().max() > 30:
                raise ValueError(
                    "GSE10072 was expected to contain RMA log2 expression values"
                )
            gene_matrix, titles = collapse_probes_to_genes(matrix, annotation)
            table = analyze_unpaired(gene_matrix, samples, titles)
        else:
            samples = classify_gse7670(metadata)
            if samples.pair.nunique() != 27:
                raise ValueError(
                    f"Expected 27 matched GSE7670 pairs; found {samples.pair.nunique()}"
                )
            matrix = matrix[samples.gsm]
            gene_matrix, titles = collapse_probes_to_genes(
                log2_transform(matrix), annotation
            )
            table = analyze_paired(gene_matrix, samples, titles)

        manifest = metadata.merge(samples, on="gsm", how="inner")
        table.to_csv(RESULTS / f"{accession}_dgea.csv", index=False)
        manifest.to_csv(RESULTS / f"{accession}_samples.csv", index=False)
        plot_volcano(table, accession)
        tables[accession] = table
        matrices[accession] = gene_matrix
        samples_by_accession[accession] = samples
        print(
            f"{accession}: {len(table):,} genes; {len(manifest)} samples; "
            f"FDR<0.05: {table.p_adj_bh.lt(0.05).sum():,}; "
            "FDR<0.05 and |log2FC|>=1: "
            f"{(table.p_adj_bh.lt(0.05) & table.log2_fold_change_tumor_vs_normal.abs().ge(1)).sum():,}"
        )

    plot_gse7670_paired_genes(
        matrices["GSE7670"],
        samples_by_accession["GSE7670"],
        tables["GSE7670"],
    )
    comparison, rho, rho_p = compare_datasets(tables)
    write_summary(tables, comparison, rho, rho_p)
    print(
        f"Cross-dataset: {comparison.significant_both.sum():,} significant in both; "
        f"{(comparison.significant_both & comparison.same_direction).sum():,} concordant"
    )
    print(f"Shared-gene Spearman rho={rho:.3f}, p={rho_p:.3g}")


if __name__ == "__main__":
    main()
