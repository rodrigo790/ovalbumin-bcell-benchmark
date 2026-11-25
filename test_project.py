import unittest

import pandas as pd

from analysis import audit, evaluate, METHODS, score_tracks
from common import DATA, RESULTS


class BCellBenchmarkTests(unittest.TestCase):
    def test_scientific_join_is_reproducible(self):
        sequence, tracks = score_tracks()
        assays = pd.read_csv(DATA / "assays.csv").fillna("")
        paired, *_ = audit(assays, tracks, sequence)
        self.assertFalse(paired[list(METHODS)].isna().any().any())
        current = evaluate(paired)
        self.assertGreaterEqual(len(current), 6)
        saved = pd.read_csv(RESULTS / "metrics.csv").fillna("")
        pd.testing.assert_frame_equal(current, saved, atol=1e-6, check_dtype=False)
