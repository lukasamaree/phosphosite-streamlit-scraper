import csv
import tempfile
import unittest
from pathlib import Path

from protein_identity_mapper import best_ncbi_doc, collect_names_from_outputs, normalize_name


class ProteinIdentityMapperTests(unittest.TestCase):
    def test_normalize_name_removes_organism_and_isoform_noise(self):
        self.assertEqual(normalize_name("Akt1 iso1 (human)"), "Akt1")
        self.assertEqual(normalize_name(" ERK2  (mouse) "), "ERK2")

    def test_collect_names_from_outputs_reads_scraper_protein_columns(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            output = Path(temp_dir) / "Akt1"
            output.mkdir()
            csv_path = output / "site.csv"
            with open(csv_path, "w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=["Protein", "Downstream protein", "Upstream protein", "Activity"],
                )
                writer.writeheader()
                writer.writerow(
                    {
                        "Protein": "Akt1",
                        "Downstream protein": "MDM2",
                        "Upstream protein": "DRD2",
                        "Activity": "ignored",
                    }
                )
                writer.writerow(
                    {
                        "Protein": "Akt1",
                        "Downstream protein": "TRAF6",
                        "Upstream protein": "DRD3",
                        "Activity": "ignored",
                    }
                )

            names = collect_names_from_outputs(temp_dir)

        self.assertEqual(names, ["Akt1", "MDM2", "DRD2", "TRAF6", "DRD3"])

    def test_best_ncbi_doc_requires_exact_symbol_or_alias_match(self):
        docs = [
            {"name": "TP53", "nomenclaturesymbol": "TP53", "otheraliases": "P53,TRP53"},
            {"name": "MDM2", "nomenclaturesymbol": "MDM2", "otheraliases": "HDM2"},
        ]

        self.assertEqual(best_ncbi_doc("MDM2", docs)["name"], "MDM2")
        self.assertEqual(best_ncbi_doc("Akt1", docs), {})


if __name__ == "__main__":
    unittest.main()
