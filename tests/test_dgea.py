import unittest

import numpy as np
import pandas as pd

from dgea import (
    analyze_paired,
    benjamini_hochberg,
    classify_gse10072,
    classify_gse7670,
    collapse_probes_to_genes,
    log2_transform,
)


class AnalysisTests(unittest.TestCase):
    def test_bh_known_values_and_order(self):
        actual = benjamini_hochberg([0.01, 0.04, 0.03, 0.002])
        np.testing.assert_allclose(actual, [0.02, 0.04, 0.04, 0.008])

    def test_gse7670_excludes_non_patient_and_incomplete_samples(self):
        metadata = pd.DataFrame(
            {
                "gsm": ["a", "b", "c", "d", "e"],
                "title": [
                    "1N: adjacent normal",
                    "1T: tumor",
                    "2N: adjacent normal",
                    "Cell line",
                    "mixture",
                ],
            }
        )
        result = classify_gse7670(metadata)
        self.assertEqual(result.pair.tolist(), ["1", "1"])

    def test_gse10072_keeps_all_tumor_and_normal_samples(self):
        metadata = pd.DataFrame(
            {
                "gsm": ["a", "b", "c", "d"],
                "title": ["Tumor_GT01", "Normal_GT01", "Tumor_GT02", "other"],
                "source": ["Adenocarcinoma", "Normal Lung Tissue", "Adenocarcinoma", "cell line"],
            }
        )
        result = classify_gse10072(metadata)
        self.assertEqual(result.gsm.tolist(), ["a", "b", "c"])
        self.assertEqual(result.group.tolist(), ["tumor", "normal", "tumor"])

    def test_log2_transform_returns_log_scale(self):
        matrix = pd.DataFrame({"s1": [1.0, 4.0], "s2": [2.0, 8.0]})
        actual = log2_transform(matrix)
        np.testing.assert_allclose(actual.to_numpy(), [[0.0, 1.0], [2.0, 3.0]])

    def test_log2_transform_rejects_nonpositive_values(self):
        with self.assertRaises(ValueError):
            log2_transform(pd.DataFrame({"s1": [0.0, 1.0]}))

    def test_probe_mapping_excludes_ambiguous_and_averages_duplicates(self):
        matrix = pd.DataFrame(
            {"s1": [2.0, 4.0, 9.0], "s2": [6.0, 8.0, 10.0]},
            index=["p1", "p2", "p3"],
        )
        annotation = pd.DataFrame(
            {
                "probe_id": ["p1", "p2", "p3"],
                "gene_symbol": ["A", "A", pd.NA],
                "gene_title": ["gene A", "gene A", "ambiguous"],
            }
        )
        genes, titles = collapse_probes_to_genes(matrix, annotation)
        self.assertEqual(genes.index.tolist(), ["A"])
        np.testing.assert_allclose(genes.loc["A"], [3.0, 7.0])
        self.assertEqual(titles.loc["A"], "gene A")

    def test_paired_analysis_reports_expected_log2_fold_change(self):
        matrix = pd.DataFrame(
            {
                "n1": [1.0, 4.0],
                "t1": [2.0, 2.0],
                "n2": [1.0, 5.0],
                "t2": [2.0, 3.0],
                "n3": [1.0, 6.0],
                "t3": [2.0, 4.0],
            },
            index=["UP", "DOWN"],
        )
        samples = pd.DataFrame(
            {
                "gsm": ["n1", "t1", "n2", "t2", "n3", "t3"],
                "pair": ["1", "1", "2", "2", "3", "3"],
                "group": ["normal", "tumor"] * 3,
            }
        )
        titles = pd.Series({"UP": "up gene", "DOWN": "down gene"})
        result = analyze_paired(matrix, samples, titles).set_index("gene_symbol")
        self.assertAlmostEqual(result.loc["UP", "log2_fold_change_tumor_vs_normal"], 1.0)
        self.assertAlmostEqual(result.loc["DOWN", "log2_fold_change_tumor_vs_normal"], -2.0)


if __name__ == "__main__":
    unittest.main()
