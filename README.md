# Paired differential gene expression in lung adenocarcinoma

For this CS 123A project, we compare lung adenocarcinoma samples with matched normal lung tissue. Comparing two tissues from the same patient helps account for baseline differences between people. We analyze each public study separately to see which expression changes recur in the second cohort.

## Results

Of the 12,502 genes tested in both studies, 3,659 have Benjamini-Hochberg adjusted p-values below 0.05 in both. Of those, 3,655 change in the same direction. Across all shared genes, estimated log2 fold changes have a Spearman rank correlation of 0.783.

Examples with higher tumor expression in both studies include `SPP1`, `MMP1`, and `COL11A1`. Examples with lower tumor expression include `AGER`, `FABP4`, and `SFTPC`. In GSE7670, `AGER` decreased in all 27 matched pairs, and `SPP1` increased in 26.

| Dataset | Matched pairs | Genes tested | BH FDR < 0.05 | FDR < 0.05 and abs(log2FC) >= 1 |
|---|---:|---:|---:|---:|
| GSE10072 | 33 | 12,502 | 6,223 | 436 |
| GSE7670 | 27 | 12,502 | 4,366 | 1,123 |

The last column counts the points highlighted in the volcano plots. A gene is significant if its BH-adjusted p-value is below 0.05. The extra cutoff of an absolute log2 fold change of 1 marks changes of at least twofold in either direction.

## Analysis design

For GSE10072, we match tumor and normal profiles that share a GEO `GT` subject ID. The primary analysis uses the 33 IDs with exactly one tumor and one normal sample. The other 41 profiles do not have a unique matched counterpart and are excluded from the paired test.

For GSE7670, we use 27 numbered tumor and adjacent-normal patient pairs. Cell lines, mixtures, commercial tissue, and other unpaired samples are excluded.

Both datasets use the Affymetrix GPL96 platform. GSE10072's processed RMA values are already on a log2 scale. We log2-transform GSE7670's positive MAS5-scaled values. We map probes using GEO's GPL96 annotation dated August 9, 2016, drop blank and ambiguous gene-symbol mappings, and average retained probes for each gene within each sample. We estimate the average tumor-minus-normal log2 expression difference and use a two-sided paired t-test for each gene. We then apply Benjamini-Hochberg correction within each study. More detail is in [`docs/data_notes.md`](docs/data_notes.md).

## Figures

- [`results/GSE10072_volcano.png`](results/GSE10072_volcano.png): effect size and adjusted significance for GSE10072
- [`results/GSE7670_volcano.png`](results/GSE7670_volcano.png): effect size and adjusted significance for GSE7670
- [`results/GSE7670_paired_genes.png`](results/GSE7670_paired_genes.png): matched-patient examples for `AGER` and `SPP1`
- [`results/cross_dataset_concordance.png`](results/cross_dataset_concordance.png): fold-change agreement between the two datasets

## Run the project

Run these commands from the repository root after cloning. The workflow was verified with Python 3.13. The first analysis run needs Internet access to download the two series matrices and GPL96 annotation from NCBI GEO.

```sh
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python scripts/run_analysis.py
python scripts/make_figures.py
python -m unittest discover -s tests
```

The analysis caches the three source files in `data/` and writes sample manifests and gene tables to `results/`. To rerun without Internet access, place the three original gzip files there first. The exact download URLs are listed in [`scripts/run_analysis.py`](scripts/run_analysis.py).

The figure script regenerates the four plots and prints the cross-study counts and Spearman correlation. A successful run returns 33 and 27 matched pairs, 6,223 and 4,366 significant genes, and the cross-study results above. The unit tests also check that the paired test uses patient IDs rather than row order.

## Repository map

- `dgea.py`: parsing, sample selection, probe collapsing, paired tests, and FDR adjustment
- `scripts/run_analysis.py`: downloads the public inputs and reproduces the result tables
- `scripts/make_figures.py`: creates the volcano, paired-patient, and concordance figures
- `tests/test_dgea.py`: checks pairing, transformations, annotation rules, statistics, and FDR adjustment
- `results/`: sample manifests, result tables, and final figures
- `presentation/`: [presentation deck](presentation/CS123A_DGEA_Presentation_Final.pptx), [shared speaker notes](presentation/speaker_notes.md), and [Veedhi's results script](presentation/veedhi_script.md)

## Presentation split

- Kareem: project motivation, datasets, sample pairing, and analysis workflow
- Veedhi: dataset results, matched-patient examples, cross-dataset validation, conclusion, and limitations
- Both: final review and questions

## Interpretation and limits

The agreement between independent cohorts supports a reproducible lung adenocarcinoma expression signature, but this is an observational class project. It uses processed microarray data from one older platform, does not adjust for clinical covariates, and does not establish that any gene causes cancer. The results are not a medical or diagnostic tool.

## Data sources

- GSE10072: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE10072
- GSE7670: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE7670
- GPL96 annotation: https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL96/annot/GPL96.annot.gz
