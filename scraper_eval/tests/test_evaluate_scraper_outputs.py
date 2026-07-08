from pathlib import Path
import unittest

from scraper_eval.evaluate_scraper_outputs import evaluate


FIXTURES = Path(__file__).resolve().parent / "fixtures"
MANIFEST = FIXTURES / "manifest.csv"


class ScraperOutputEvaluationTests(unittest.TestCase):
    def test_correct_outputs_score_100(self):
        summary = evaluate(MANIFEST, FIXTURES / "correct_ids.csv", FIXTURES / "correct_outputs")
        result = summary["results"][0]
        self.assertEqual(summary["status"], "passed")
        self.assertEqual(summary["wrong_protein_rate"], 0)
        self.assertEqual(result["identity_score"], 1)
        self.assertEqual(result["ptm_score"], 1)
        self.assertEqual(result["final_score"], 100)

    def test_wrong_protein_id_is_critical_and_zeroes_score(self):
        summary = evaluate(MANIFEST, FIXTURES / "wrong_id_ids.csv", FIXTURES / "wrong_outputs")
        result = summary["results"][0]
        self.assertEqual(summary["status"], "failed")
        self.assertEqual(summary["wrong_protein_rate"], 1)
        self.assertIs(result["critical_identity_failure"], True)
        self.assertEqual(result["identity_score"], 0)
        self.assertEqual(result["ptm_score"], 0)
        self.assertEqual(result["final_score"], 0)

    def test_wrong_organism_is_critical(self):
        summary = evaluate(MANIFEST, FIXTURES / "wrong_organism_ids.csv", FIXTURES / "wrong_outputs")
        result = summary["results"][0]
        self.assertIs(result["critical_identity_failure"], True)
        self.assertEqual(result["final_score"], 0)
        self.assertTrue(any("wrong organism" in error for error in result["errors"]))

    def test_wrong_uniprot_is_critical(self):
        summary = evaluate(MANIFEST, FIXTURES / "wrong_uniprot_ids.csv", FIXTURES / "wrong_outputs")
        result = summary["results"][0]
        self.assertIs(result["critical_identity_failure"], True)
        self.assertEqual(result["final_score"], 0)
        self.assertTrue(any("wrong UniProt" in error for error in result["errors"]))

    def test_duplicate_ptm_rows_fail_ptm_only(self):
        summary = evaluate(MANIFEST, FIXTURES / "correct_ids.csv", FIXTURES / "duplicate_outputs")
        result = summary["results"][0]
        self.assertEqual(result["identity_score"], 1)
        self.assertLess(result["ptm_score"], 1)
        self.assertLess(result["final_score"], 100)
        self.assertEqual(result["duplicate_row_count"], 1)

    def test_schema_errors_fail_ptm_only(self):
        summary = evaluate(MANIFEST, FIXTURES / "correct_ids.csv", FIXTURES / "schema_bad_outputs")
        result = summary["results"][0]
        self.assertEqual(result["identity_score"], 1)
        self.assertLess(result["ptm_score"], 1)
        self.assertTrue(result["schema_errors"])


if __name__ == "__main__":
    unittest.main()
