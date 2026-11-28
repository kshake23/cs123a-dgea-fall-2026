# How we picked the samples

## GSE10072

GEO reports 107 final expression profiles: 58 lung tumors and 49 normal lung samples. The titles include subject IDs such as `GT00006`, and many IDs have both a tumor and normal sample. We label tissue from `source_name_ch1`, then pair only IDs that have exactly one sample of each type. The analysis uses the 33 complete pairs (66 samples) and leaves out 41 profiles without a unique matched counterpart. GEO also notes that duplicate or triplicate QC samples were averaged while preparing the final expression values. We use the processed RMA log2 values as provided.

## GSE7670

GEO has 66 samples total. Twenty-seven patients have tumor and adjacent-normal samples identified by matching numeric labels, such as `1N` and `1T`. We keep those 54 samples and leave out the other 12 samples, which include cell lines, mixtures, and commercial normal tissue. The series matrix contains positive MAS5-scaled values, which we convert to log2 before averaging probes and testing.

Both datasets use the Affymetrix Human Genome U133A platform, GPL96, and remain separate during testing. The workflow uses GEO's processed series-matrix values rather than raw CEL files.

## GPL96 annotation and probe handling

The gene mapping comes from GEO's `GPL96.annot.gz` file, whose header reports `!Annotation_date = Aug 09 2016`. A probe is retained only when its gene-symbol field is nonblank, is not `---`, and contains one symbol rather than a `///`-separated ambiguous mapping. After filtering, expression values for multiple probes with the same gene symbol are averaged for each sample. The paired statistical tests and Benjamini-Hochberg correction are then run at the gene level.

## Reproduction inputs

The three gzip inputs downloaded from NCBI GEO for the verified run had these SHA-256 checksums:

```text
88e0b22362bac779eb220b3b185c80faa6510a92b9358eaad159a561ab4351c4  GPL96.annot.gz
68fc224d6c23b1802ddc2223a0433c8cf25f11cdb9606e33736f4a0c098496af  GSE10072_series_matrix.txt.gz
f321f33535c8b929f08499e12e3c2a950030b84bfeb47de8ba8bbd32b934dad7  GSE7670_series_matrix.txt.gz
```
