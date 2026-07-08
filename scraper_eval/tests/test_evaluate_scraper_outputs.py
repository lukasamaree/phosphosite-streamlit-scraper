from pathlib import Path
import tempfile
import unittest

from scraper_eval.evaluate_scraper_outputs import evaluate
from validate_scrape_outputs import build_summary


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
        self.assertEqual(result["failure_class"], "VALID_SCRAPE")
        self.assertEqual(result["recommended_action"], "use_output")
        self.assertEqual(summary["classification"]["usable_scrape_rate"], 1)

    def test_wrong_protein_id_is_critical_and_zeroes_score(self):
        summary = evaluate(MANIFEST, FIXTURES / "wrong_id_ids.csv", FIXTURES / "wrong_outputs")
        result = summary["results"][0]
        self.assertEqual(summary["status"], "failed")
        self.assertEqual(summary["wrong_protein_rate"], 1)
        self.assertIs(result["critical_identity_failure"], True)
        self.assertEqual(result["identity_score"], 0)
        self.assertEqual(result["ptm_score"], 0)
        self.assertEqual(result["final_score"], 0)
        self.assertEqual(result["failure_class"], "IDENTITY_FAILURE")
        self.assertEqual(result["recommended_action"], "stop_and_resolve_identity")

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
        self.assertEqual(result["failure_class"], "SCHEMA_FAILURE")

    def test_placeholder_rows_fail_ptm_only(self):
        summary = evaluate(MANIFEST, FIXTURES / "correct_ids.csv", FIXTURES / "placeholder_outputs")
        result = summary["results"][0]
        self.assertEqual(result["identity_score"], 1)
        self.assertLess(result["ptm_score"], 1)
        self.assertEqual(result["placeholder_row_count"], 1)
        self.assertEqual(result["failure_class"], "CLOUDFLARE_LIKELY")
        self.assertIn("placeholder_rows_detected", result["failure_signals"])
        self.assertEqual(summary["classification"]["cloudflare_likely_rate"], 1)
        self.assertTrue(any("placeholder PTM rows" in error for error in result["errors"]))

    def test_validator_uses_curated_ids_under_requested_root(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            curated_dir = root / "curated_protein_ids"
            curated_dir.mkdir()
            (curated_dir / "resolved_protein_ids.csv").write_text(
                "protein_name,protein_id,url,organism,source\n"
                "AKT1,570,https://www.phosphosite.org/proteinAction.action?id=570&showAllSites=true,human,fixture\n",
                encoding="utf-8",
            )

            summary = build_summary(root)

        self.assertTrue(summary["curated_ids"]["exists"])
        self.assertEqual(summary["curated_ids"]["rows"], 1)
        self.assertEqual(summary["curated_ids"]["path"], "curated_protein_ids/resolved_protein_ids.csv")


if __name__ == "__main__":
    unittest.main()
