# How we picked the samples

## GSE10072

GEO reports 107 final expression profiles: 58 lung tumors and 49 normal lung samples. The titles include subject IDs such as `GT00006`, and many IDs have both a tumor and normal sample. We label tissue from `source_name_ch1`, then pair only IDs that have exactly one sample of each type. The current workflow finds 33 complete pairs (66 samples) and leaves out 41 samples with no unique complete match. We use a two-sided paired t-test on those pairs. This is a conservative paired subset; it does not use the unmatched samples or claim that all 107 profiles can be paired. GEO also notes that duplicate or triplicate QC samples were averaged as part of preparing the final expression values.

## GSE7670

GEO has 66 samples total. Twenty-seven patients have tumor and adjacent-normal samples identified by matching numeric labels, such as `1N` and `1T`. We keep those 54 samples and leave out the other 12 samples, which include cell lines, mixtures, and commercial normal tissue. The analysis uses a two-sided paired t-test for each probe.

Both datasets use GPL96, but they stay separate in the analysis. We use GEO's processed series-matrix values rather than raw CEL files. Results are for probe sets; gene-symbol mapping and combining multiple probes are not included yet.
