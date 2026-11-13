# Lung cancer gene expression project (work in progress)

We’re building a small Python tool to compare gene activity in lung tumor and normal tissue for our CS 123A project. It currently runs on two public GEO datasets. We keep them separate because they come from different studies. The code and first results are ready to review, but we haven’t decided what to present yet.

## Current status

- Kareem’s analysis code runs on GSE10072 and GSE7670.
- GSE10072 currently uses 33 complete tumor/normal pairs; GSE7670 uses 27 pairs.
- The output tables and sample lists are in `results/`; some are checked into the repo.
- Veedhi’s review of GSE7670 and presentation figures are still pending.
- The team still needs to review the findings, decide what to show, and finish the presentation.

## Get started

```sh
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python scripts/run_analysis.py
python -m unittest discover -s tests
```

The script downloads each dataset’s GEO expression file into `data/`, then saves the sample lists, test results, and plots in `results/`. GEO has already processed the expression values, and we analyze each study separately. Plots are created on your computer; we’ll decide which ones to use after checking the results.

## What the analysis does

For GSE10072, the current version matches tumor and normal samples that share a GEO `GT` ID, then compares the 33 complete pairs. The 41 other samples don’t have a unique match, so they’re left out of this test. For GSE7670, it uses 27 numbered tumor/normal pairs and leaves out cell lines, mixtures, and other unpaired samples. The code tests each array probe and adjusts for the many tests it runs. We still need to check these choices and results as a team.

The current tables show results for array probes. Several probes can point to the same gene, and we haven’t decided how to combine them yet. So these aren’t our final gene-level results. We should check the sample lists and tables together before using results in the presentation.

## Data sources

- GSE10072: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE10072
- GSE7670: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE7670

Both use the Affymetrix GPL96 platform. More details about which samples are included are in `docs/data_notes.md`.

This is a class project for learning, not a medical tool. We should write final project prose in our own words and follow the course rules on AI use.
