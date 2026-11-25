# Paired differential gene expression in lung adenocarcinoma

> **Working draft (November 24, 2026):** The analysis and current figures are in the repo, but we still need to go over the results and finish the written explanation together. We may update the interpretation and presentation after that review.

For this CS 123A project, we compare tumor and normal lung samples from two public studies. We analyze each study on its own, then check which expression changes show up in both.

## What the current results show

In the current run, 3,659 genes are significant in both studies using a Benjamini-Hochberg false discovery rate (FDR) below 0.05. Of those, 3,655 change in the same direction. The estimated fold changes across the two studies have a Spearman correlation of 0.783.

Some examples with higher tumor expression in both studies are `SPP1`, `MMP1`, and `COL11A1`. Some examples with lower tumor expression are `AGER`, `FABP4`, and `SFTPC`. We still need to review these results and decide how to describe them in the presentation.

## Current results for review

| Dataset | Matched pairs | Genes tested | BH FDR < 0.05 | FDR < 0.05 and abs(log2FC) >= 1 |
|---|---:|---:|---:|---:|
| GSE10072 | 33 | 12,502 | 6,223 | 436 |
| GSE7670 | 27 | 12,502 | 4,366 | 1,123 |

The last column only marks larger changes in the volcano plots. For the analysis, a gene is considered significant when its BH FDR is below 0.05. These counts come from the current code and should be reviewed together before we use them as final presentation claims.

## Analysis design

For GSE10072, we match tumor and normal profiles that share a GEO `GT` subject ID. The primary analysis uses the 33 IDs with exactly one tumor and one normal sample. The other 41 profiles do not have a unique matched counterpart and are excluded from the paired test.

For GSE7670, we use 27 numbered tumor and adjacent-normal patient pairs. Cell lines, mixtures, commercial tissue, and other unpaired samples are excluded.

Both datasets use the Affymetrix GPL96 platform. Probes are mapped using GEO's GPL96 annotation dated August 9, 2016. Blank and ambiguous gene-symbol mappings are removed, and multiple retained probes for the same gene are averaged. Each gene is tested with a two-sided paired t-test, followed by Benjamini-Hochberg correction across genes. More detail is in [`docs/data_notes.md`](docs/data_notes.md).

## Figures currently in the repo

- [`results/GSE10072_volcano.png`](results/GSE10072_volcano.png): effect size and adjusted significance for GSE10072
- [`results/GSE7670_volcano.png`](results/GSE7670_volcano.png): effect size and adjusted significance for GSE7670
- [`results/GSE7670_paired_genes.png`](results/GSE7670_paired_genes.png): matched-patient examples for `AGER` and `SPP1`
- [`results/cross_dataset_concordance.png`](results/cross_dataset_concordance.png): fold-change agreement between the two datasets

## Run the project

```sh
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python scripts/run_analysis.py
python scripts/make_figures.py
python -m unittest discover -s tests
```

The analysis script downloads the GEO series matrices and GPL96 annotation into `data/`, then writes sample manifests and gene-level results to `results/`. The figure script reads those outputs and creates all four presentation figures.

## Repository map

- `dgea.py`: parsing, sample selection, probe collapsing, paired tests, and FDR adjustment
- `scripts/run_analysis.py`: downloads the public inputs and reproduces the result tables
- `scripts/make_figures.py`: creates the volcano, paired-patient, and concordance figures
- `tests/test_dgea.py`: checks pairing, transformations, annotation rules, statistics, and FDR adjustment
- `results/`: sample manifests, result tables, and final figures
- `presentation/`: presentation deck and speaker notes

## Presentation split

- Kareem: project motivation, datasets, sample pairing, and analysis workflow
- Veedhi: dataset results, matched-patient examples, cross-dataset validation, conclusion, and limitations
- Both: final review and questions

## Still to finish

- Review the results and figures together and agree on what to highlight.
- Finish the written explanation and update this README if our review changes anything.
- Each of us needs to submit the GitHub repository link before presenting.
- Veedhi's latest email says Canvas showed our DGEA presentation slot on Monday, November 30. Check Canvas again before we present.

## Interpretation and limits

The agreement between independent cohorts supports a reproducible lung adenocarcinoma expression signature, but this is an observational class project. It uses processed microarray data from one older platform, does not adjust for clinical covariates, and does not establish that any gene causes cancer. The results are not a medical or diagnostic tool.

## Data sources

- GSE10072: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE10072
- GSE7670: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE7670
- GPL96 annotation: https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL96/annot/GPL96.annot.gz
