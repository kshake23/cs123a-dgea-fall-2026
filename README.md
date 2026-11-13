# Lung cancer gene expression project

This is our CS 123A project for comparing lung tumor and normal tissue expression with Python. We use two public GEO datasets and analyze each one on its own because the samples come from different studies.

## Get started

```sh
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python scripts/run_analysis.py
python -m unittest discover -s tests
```

The script downloads each dataset's GEO series matrix into `data/`, then saves the sample lists, test results, and volcano plots in `results/`. The data are already processed expression values from GEO, so the script does not treat them as raw measurements or combine the studies.

## What the analysis does

For GSE10072, it pairs tumor and normal samples with the same GEO `GT` identifier and uses a two-sided paired t-test. The current run keeps 33 complete pairs and leaves out 41 samples without a unique matching pair. For GSE7670, it keeps the 27 numbered tumor/normal patient pairs and also uses a paired t-test. Samples such as cell lines and mixtures are left out of that paired comparison. The code adjusts p-values across probe sets with the Benjamini-Hochberg method.

The results are at the microarray probe level. A probe can map to a gene, and more than one probe can map to the same gene, so we leave symbol-level combining for a later step with an agreed rule. Check the sample lists and results before using them in the presentation.

## Data sources

- GSE10072: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE10072
- GSE7670: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE7670

Both use the Affymetrix GPL96 platform. More details about which samples are included are in `docs/data_notes.md`.

This is a class project for learning, not a medical tool. We should write final project prose in our own words and follow the course rules on AI use.
