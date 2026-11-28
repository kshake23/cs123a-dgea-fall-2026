# CS 123A presentation notes

## Slide 1: Paired differential gene expression in lung adenocarcinoma

**Speaker: Kareem, about 30 seconds**

We looked at two public lung adenocarcinoma datasets to ask a simple question: which expression changes are reproducible when tumor tissue is compared with matched normal lung from the same patient? The main point is not just finding significant genes in one dataset. It is checking whether the direction and size of those changes hold up in a second study.

## Slide 2: Why matched samples matter

**Speaker: Kareem, about 1 minute**

Gene expression varies a lot between people. If tumor and normal samples come from different people, some of the difference may be ordinary person-to-person variation. A paired design compares each patient's tumor with their own normal tissue, so the test focuses on the within-patient change. Using two independent cohorts then gives us a separate check on reproducibility.

## Slide 3: Two independent cohorts, one shared platform

**Speaker: Kareem, about 1 minute**

GSE10072 has 33 complete tumor-normal pairs after we require exactly one sample of each type for a subject ID. GSE7670 has 27 numbered tumor and adjacent-normal pairs after excluding cell lines, mixtures, commercial tissue, and unpaired samples. Both studies use the Affymetrix GPL96 platform, which lets us apply the same gene-mapping rules while still analyzing the studies separately.

## Slide 4: From GEO files to gene-level tests

**Speaker: Kareem, about 2 minutes**

The workflow starts from GEO's processed series matrices and the GPL96 annotation. We keep only unambiguous gene-symbol mappings, average multiple retained probes for the same gene, and build the matched sample pairs. For each gene we run a two-sided paired t-test and then apply Benjamini-Hochberg correction across all genes. We use FDR below 0.05 for statistical significance. The extra absolute log2 fold-change cutoff of 1 is only used to highlight larger effects in the volcano plots.

## Slide 5: GSE10072 shows widespread tumor-associated change

**Speaker: Veedhi, about 1 minute**

In GSE10072 we tested 12,502 genes across 33 matched pairs. A total of 6,223 genes have adjusted p-values below 0.05. The orange points are the 436 genes that also have an absolute log2 fold change of at least 1. Examples with higher tumor expression include SPP1, MMP1, and COL11A1. Examples with lower tumor expression include AGER, FABP4, and SFTPC.

## Slide 6: GSE7670: second cohort

**Speaker: Veedhi, about 1 minute**

In the second study, 4,366 of the same 12,502 genes are significant across 27 matched pairs. There are 1,123 orange points with both FDR below 0.05 and an absolute log2 fold change of at least 1. The exact counts differ across studies, which is expected, but many of the strongest genes and their directions are the same.

## Slide 7: The paired view makes the result tangible

**Speaker: Veedhi, about 1 minute**

These plots show the individual GSE7670 patients rather than only summary statistics. Each thin line connects one patient's normal and tumor samples. AGER drops in nearly every pair, while SPP1 rises in nearly every pair. The thick line connects the group means. This is a useful check that the result is not being driven by only one or two extreme samples.

## Slide 8: The two studies agree gene by gene

**Speaker: Veedhi, about 1.5 minutes**

Here each point is one gene. Its horizontal position is the GSE10072 fold change and its vertical position is the GSE7670 fold change. The overall Spearman correlation is 0.783. There are 3,659 genes that are significant in both studies, and 3,655 of them change in the same direction. The labeled genes show examples at both ends of the plot. This agreement is stronger evidence than a single-dataset result.

## Slide 9: What we can conclude

**Speaker: Veedhi, about 1 minute, then both for questions**

The two independent cohorts support a reproducible lung adenocarcinoma expression signature. The strongest repeated changes include higher SPP1, MMP1, and COL11A1 and lower AGER, FABP4, and SFTPC in tumors. The limits are important: these are processed microarray data from an older platform, we did not adjust for clinical covariates, and association does not establish causation. The full code, sample manifests, results, tests, and figures are in the repository so the analysis can be rerun.

## Timing and handoff

- Kareem: slides 1 through 4. Rehearse to the 5–7 minute window.
- Veedhi: slides 5 through 9. [Timed script and Q&A](veedhi_script.md) target about 6 minutes.
- Both: questions

Practice the handoff after slide 4: "Veedhi will walk through what we found in each dataset and how well the results reproduced."
