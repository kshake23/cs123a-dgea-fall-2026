# How we picked the samples

## GSE10072

GEO reports 107 final expression profiles: 58 lung tumors and 49 normal lung samples. The titles include subject IDs such as `GT00006`, and many IDs have both a tumor and normal sample. We label tissue from `source_name_ch1`, then pair only IDs that have exactly one sample of each type. The primary analysis uses the 33 complete pairs (66 samples) and leaves out 41 profiles without a unique matched counterpart. This keeps the paired design described in the project plan. An all-sample Welch test is not part of this first increment and would be treated as a secondary sensitivity analysis if added later. GEO also notes that duplicate or triplicate QC samples were averaged while preparing the final expression values.

## GSE7670

GEO has 66 samples total. Twenty-seven patients have tumor and adjacent-normal samples identified by matching numeric labels, such as `1N` and `1T`. We keep those 54 samples and leave out the other 12 samples, which include cell lines, mixtures, and commercial normal tissue. After probe-to-gene mapping, the analysis uses a two-sided paired t-test for each gene.

Both datasets use the Affymetrix Human Genome U133A platform, GPL96, and remain separate during testing. The workflow uses GEO's processed series-matrix values rather than raw CEL files.

## GPL96 annotation and probe handling

The gene mapping comes from GEO's `GPL96.annot.gz` file, whose header reports `!Annotation_date = Aug 09 2016`. A probe is retained only when its gene-symbol field is nonblank, is not `---`, and contains one symbol rather than a `///`-separated ambiguous mapping. After filtering, expression values for multiple probes with the same gene symbol are averaged for each sample. The paired statistical tests and Benjamini-Hochberg correction are then run at the gene level.
