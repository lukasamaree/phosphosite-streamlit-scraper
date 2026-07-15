import csv
import tempfile
import unittest
from pathlib import Path

from protein_identity_mapper import (
    best_ncbi_doc,
    collect_names_from_outputs,
    enrich_outputs,
    is_probable_protein_token,
    normalize_name,
)


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

    def test_rejects_relationship_text_as_protein_token(self):
        self.assertFalse(is_probable_protein_token("Effects of modification on Akt1"))
        self.assertFalse(is_probable_protein_token("Regulatory protein"))
        self.assertTrue(is_probable_protein_token("Akt1"))
        self.assertTrue(is_probable_protein_token("KPNB1"))

    def test_enrich_outputs_adds_gene_and_uniprot_columns(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            output_dir = root / "Akt1"
            output_dir.mkdir()
            input_csv = output_dir / "site.csv"
            lookup_csv = root / "protein_identity_lookup.csv"
            enriched_root = root / "enriched"

            with open(input_csv, "w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=["Protein", "Downstream protein", "Upstream protein"],
                )
                writer.writeheader()
                writer.writerow({"Protein": "Akt1", "Downstream protein": "MDM2", "Upstream protein": "DRD2"})

            with open(lookup_csv, "w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=[
                        "raw_name",
                        "normalized_query",
                        "canonical_gene",
                        "uniprot_accession",
                        "recommended_name",
                        "organism",
                        "aliases",
                        "confidence",
                        "match_status",
                        "sources",
                        "hgnc_id",
                        "ncbi_gene_id",
                        "ensembl_gene_id",
                    ],
                )
                writer.writeheader()
                writer.writerow({"raw_name": "Akt1", "canonical_gene": "AKT1", "uniprot_accession": "P31749", "confidence": "high", "sources": "UniProt", "aliases": "Akt1;PKB"})
                writer.writerow({"raw_name": "MDM2", "canonical_gene": "MDM2", "uniprot_accession": "Q00987", "confidence": "high", "sources": "UniProt", "aliases": "MDM2"})
                writer.writerow({"raw_name": "DRD2", "canonical_gene": "DRD2", "uniprot_accession": "P14416", "confidence": "high", "sources": "UniProt", "aliases": "DRD2"})

            summaries = enrich_outputs(root, lookup_csv, enriched_root)
            enriched_csv = enriched_root / "Akt1" / "site.csv"

            self.assertEqual(len(summaries), 1)
            with open(enriched_csv, "r", encoding="utf-8", newline="") as handle:
                row = next(csv.DictReader(handle))

        self.assertEqual(row["protein_canonical_gene"], "AKT1")
        self.assertEqual(row["protein_uniprot_accession"], "P31749")
        self.assertEqual(row["downstream_protein_canonical_gene"], "MDM2")
        self.assertEqual(row["downstream_protein_uniprot_accession"], "Q00987")
        self.assertEqual(row["upstream_protein_canonical_gene"], "DRD2")
        self.assertEqual(row["upstream_protein_uniprot_accession"], "P14416")

    def test_enrich_outputs_leaves_blank_or_non_protein_cells_blank(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            output_dir = root / "Akt1"
            output_dir.mkdir()
            input_csv = output_dir / "site.csv"
            lookup_csv = root / "protein_identity_lookup.csv"
            enriched_root = root / "enriched"

            with open(input_csv, "w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=["Protein", "Downstream protein", "Upstream protein"],
                )
                writer.writeheader()
                writer.writerow(
                    {
                        "Protein": "Akt1",
                        "Downstream protein": "Effects of modification on Akt1",
                        "Upstream protein": "",
                    }
                )

            with open(lookup_csv, "w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=[
                        "raw_name",
                        "normalized_query",
                        "canonical_gene",
                        "uniprot_accession",
                        "recommended_name",
                        "organism",
                        "aliases",
                        "confidence",
                        "match_status",
                        "sources",
                        "hgnc_id",
                        "ncbi_gene_id",
                        "ensembl_gene_id",
                    ],
                )
                writer.writeheader()
                writer.writerow({"raw_name": "Akt1", "canonical_gene": "AKT1", "uniprot_accession": "P31749", "confidence": "high", "sources": "UniProt", "aliases": "Akt1"})
                writer.writerow({"raw_name": "SCN11A", "canonical_gene": "SCN11A", "uniprot_accession": "Q9UKU5", "confidence": "high", "sources": "UniProt", "aliases": "Effects of modification on Akt1"})

            enrich_outputs(root, lookup_csv, enriched_root)
            enriched_csv = enriched_root / "Akt1" / "site.csv"
            with open(enriched_csv, "r", encoding="utf-8", newline="") as handle:
                row = next(csv.DictReader(handle))

        self.assertEqual(row["protein_canonical_gene"], "AKT1")
        self.assertEqual(row["downstream_protein_canonical_gene"], "")
        self.assertEqual(row["downstream_protein_uniprot_accession"], "")
        self.assertEqual(row["upstream_protein_canonical_gene"], "")
        self.assertEqual(row["upstream_protein_uniprot_accession"], "")

    def test_enrich_outputs_strips_stale_identity_columns_before_recomputing(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            output_dir = root / "Akt1"
            output_dir.mkdir()
            input_csv = output_dir / "site.csv"
            lookup_csv = root / "protein_identity_lookup.csv"
            enriched_root = root / "enriched"

            with open(input_csv, "w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=[
                        "Protein",
                        "Downstream protein",
                        "Upstream protein",
                        "downstream_protein_canonical_gene",
                        "downstream_protein_uniprot_accession",
                    ],
                )
                writer.writeheader()
                writer.writerow(
                    {
                        "Protein": "Akt1",
                        "Downstream protein": "",
                        "Upstream protein": "",
                        "downstream_protein_canonical_gene": "SCN11A",
                        "downstream_protein_uniprot_accession": "Q9UKU5",
                    }
                )

            with open(lookup_csv, "w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=[
                        "raw_name",
                        "normalized_query",
                        "canonical_gene",
                        "uniprot_accession",
                        "recommended_name",
                        "organism",
                        "aliases",
                        "confidence",
                        "match_status",
                        "sources",
                        "hgnc_id",
                        "ncbi_gene_id",
                        "ensembl_gene_id",
                    ],
                )
                writer.writeheader()
                writer.writerow({"raw_name": "Akt1", "canonical_gene": "AKT1", "uniprot_accession": "P31749", "confidence": "high", "sources": "UniProt", "aliases": "Akt1"})
                writer.writerow({"raw_name": "SCN11A", "canonical_gene": "SCN11A", "uniprot_accession": "Q9UKU5", "confidence": "high", "sources": "UniProt", "aliases": "SCN11A"})

            enrich_outputs(root, lookup_csv, enriched_root)
            enriched_csv = enriched_root / "Akt1" / "site.csv"
            with open(enriched_csv, "r", encoding="utf-8", newline="") as handle:
                row = next(csv.DictReader(handle))

        self.assertEqual(row["protein_canonical_gene"], "AKT1")
        self.assertEqual(row["downstream_protein_canonical_gene"], "")
        self.assertEqual(row["downstream_protein_uniprot_accession"], "")
        self.assertEqual(row["upstream_protein_canonical_gene"], "")
        self.assertEqual(row["upstream_protein_uniprot_accession"], "")


if __name__ == "__main__":
    unittest.main()
