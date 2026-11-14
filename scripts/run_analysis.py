from pathlib import Path
import gzip
import sys
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from dgea import (
    analyze_paired,
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
ANNOTATION_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL96/annot/GPL96.annot.gz"
)


def download_gzip(path, url, required_text):
    if not path.exists():
        urllib.request.urlretrieve(url, path)
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        header = "".join(handle.readline() for _ in range(35))
    if required_text not in header:
        raise ValueError(f"{path} failed its source validation check")
    return path


def main():
    DATA.mkdir(exist_ok=True)
    RESULTS.mkdir(exist_ok=True)

    annotation_path = download_gzip(
        DATA / "GPL96.annot.gz",
        ANNOTATION_URL,
        "!Annotation_platform = GPL96",
    )
    annotation = read_platform_annotation(annotation_path)

    for accession, url in SOURCES.items():
        matrix_path = download_gzip(
            DATA / f"{accession}_series_matrix.txt.gz",
            url,
            f'!Series_geo_accession\t"{accession}"',
        )
        matrix, metadata = read_series_matrix(matrix_path)

        if accession == "GSE10072":
            samples = classify_gse10072(metadata)
            if samples.pair.nunique() != 33:
                raise ValueError(
                    f"Expected 33 matched GSE10072 pairs; found {samples.pair.nunique()}"
                )
            matrix = matrix[samples.gsm]
            if matrix.max().max() > 30:
                raise ValueError(
                    "GSE10072 was expected to contain RMA log2 expression values"
                )
        else:
            samples = classify_gse7670(metadata)
            if samples.pair.nunique() != 27:
                raise ValueError(
                    f"Expected 27 matched GSE7670 pairs; found {samples.pair.nunique()}"
                )
            matrix = log2_transform(matrix[samples.gsm])

        gene_matrix, gene_titles = collapse_probes_to_genes(matrix, annotation)
        table = analyze_paired(gene_matrix, samples, gene_titles)
        manifest = metadata.merge(samples, on="gsm", how="inner")

        table.to_csv(RESULTS / f"{accession}_dgea.csv", index=False)
        manifest.to_csv(RESULTS / f"{accession}_samples.csv", index=False)
        print(
            f"{accession}: {samples.pair.nunique()} matched pairs; "
            f"{len(table):,} genes tested; "
            f"FDR < 0.05: {table.p_adj_bh.lt(0.05).sum():,}"
        )


if __name__ == "__main__":
    main()
