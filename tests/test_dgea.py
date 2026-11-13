import unittest
import numpy as np
import pandas as pd

from dgea import benjamini_hochberg, classify_gse7670, classify_gse10072


class AnalysisTests(unittest.TestCase):
    def test_bh_known_values_and_order(self):
        actual = benjamini_hochberg([.01, .04, .03, .002])
        np.testing.assert_allclose(actual, [.02, .04, .04, .008])

    def test_gse7670_excludes_non_patient_samples(self):
        metadata = pd.DataFrame({"gsm": ["a", "b", "c", "d", "e"], "title": [
            "1N: adjacent normal", "1T: tumor", "2N: adjacent normal", "Cell line", "mixture"]})
        result = classify_gse7670(metadata)
        self.assertEqual(result.pair.tolist(), ["1", "1"])

    def test_gse10072_keeps_only_complete_id_matched_pairs(self):
        metadata = pd.DataFrame({"gsm": ["a", "b", "c"],
            "title": ["Lung Tumor_GT01", "Normal Lung_GT01", "Lung Tumor_GT02"],
            "source": ["Adenocarcinoma", "Normal Lung Tissue", "Adenocarcinoma"]})
        result = classify_gse10072(metadata)
        self.assertEqual(result.pair.unique().tolist(), ["GT01"])
        self.assertEqual(set(result.group), {"normal", "tumor"})


if __name__ == "__main__":
    unittest.main()
