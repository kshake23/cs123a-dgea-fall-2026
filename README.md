# Lung cancer gene expression project

We’re building a small Python tool to compare gene activity in lung tumor and normal tissue for our CS 123A project. It currently runs on two public GEO datasets. We keep them separate because they come from different studies. The code and first results are ready to review, but we haven’t decided what to present yet.

## Current status

- Kareem’s analysis code runs on GSE10072 and GSE7670.
- The primary methods use paired tests with 33 complete GSE10072 pairs and 27 complete GSE7670 pairs.
- Probe sets are mapped to gene symbols before testing.
- Regenerated results contain 12,502 tested genes per dataset.
- The team still needs to review the findings and finish the presentation.

## Get started

```sh
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python scripts/run_analysis.py
python scripts/make_figures.py
python -m unittest discover -s tests
```

The scripts download each dataset’s GEO expression file and the GPL96 annotation into `data/`, then save sample lists, gene-level test results, and figures in `results/`.

## What the analysis does

For GSE10072, the primary analysis matches tumor and normal samples that share a GEO `GT` ID, then compares the 33 complete pairs. The other 41 profiles do not have a unique matched counterpart, so they are excluded from the paired test. An all-sample Welch test could be added later as a secondary sensitivity analysis if the team decides it is useful. GSE7670 uses 27 numbered tumor/adjacent-normal pairs and excludes cell lines, mixtures, and other unpaired samples.

The workflow uses GEO’s GPL96 annotation dated August 9, 2016. Blank mappings, `---`, and mappings containing multiple symbols separated by `///` are excluded. When multiple retained probes map to the same gene symbol, their expression values are averaged before the paired test. Benjamini-Hochberg adjustment is applied to the gene-level p-values.

## Reproduced results

The full workflow was rerun twice from the same public inputs. Both runs produced identical output files.

| Dataset | Matched pairs | Genes tested | BH FDR < 0.05 |
|---|---:|---:|---:|
| GSE10072 | 33 | 12,502 | 6,223 |
| GSE7670 | 27 | 12,502 | 4,366 |

## Data sources

- GSE10072: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE10072
- GSE7670: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE7670
- GPL96 annotation: https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL96/annot/GPL96.annot.gz

Both use the Affymetrix GPL96 platform. More details about which samples are included are in `docs/data_notes.md`.

This is a class project for learning, not a medical or diagnostic tool.
