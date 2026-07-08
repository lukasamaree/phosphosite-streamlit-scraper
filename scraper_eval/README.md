# Scraper Output Evaluation

This framework evaluates the scraper's actual outputs, not agent prompt quality.

It is identity-first:

1. Compare resolved protein identity against an expected manifest.
2. Treat wrong protein ID, organism, or UniProt/SwissProt accession as critical failures.
3. Skip PTM scoring when identity fails.
4. Score outputs as:

```text
final_score = identity_score * ptm_score * 100
```

The summary also reports:

```text
wrong_protein_rate = critical_identity_failures / evaluated_proteins
```

## Manifest

Expected manifest columns:

```csv
query,expected_protein_id,expected_organism,expected_uniprot,output_dir,required_sites
```

`required_sites` is a semicolon-separated list, for example:

```text
Ser15;Ser20
```

## Run

```bash
.venv312/bin/python -m scraper_eval.evaluate_scraper_outputs \
  --manifest scraper_eval/tests/fixtures/manifest.csv \
  --curated-ids scraper_eval/tests/fixtures/correct_ids.csv \
  --output-root scraper_eval/tests/fixtures/correct_outputs
```

The evaluator prints a final `EVAL_JSON: ...` line so another agent can parse results.

## Tests

```bash
.venv312/bin/python -m unittest discover scraper_eval/tests
```
