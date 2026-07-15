# PhosphoSitePlus Streamlit Scraper

Streamlit dashboard and reusable command-line tools for resolving human PhosphoSitePlus protein URL IDs and scraping phosphorylation/site data.

## What It Does

- Accepts protein symbols such as `AKT1`, `TP53`, `EGFR`, and `MAPK1`.
- Searches PhosphoSitePlus and selects the human protein page.
- Reads the protein ID from the final `proteinAction.action?id=...` URL.
- Saves curated protein IDs in `curated_protein_ids/`.
- Reuses saved IDs so known proteins do not need to be looked up again.
- Scrapes resolved protein IDs with conservative randomized delays and live logs.

## Streamlit App

Run locally:

```bash
.venv312/bin/python -m streamlit run streamlit_phospho_scraper.py
```

Or use:

```bash
./app.sh
```

Deploy on Streamlit Community Cloud with:

```text
Main file path: streamlit_phospho_scraper.py
```

## Curated ID Cache

The app and tool wrapper write durable ID curation files here:

```text
curated_protein_ids/resolved_protein_ids.csv
curated_protein_ids/lookup_state.json
curated_protein_ids/protein_ids.txt
```

These files are intentionally tracked so curated protein IDs can travel with the repo. Bulky scrape outputs and generated CSV files are ignored.

## Agent Tool

Resolve and save missing human protein IDs:

```bash
.venv312/bin/python phosphosite_agent_tool.py resolve --protein-names-file protein_names.txt
```

Resolve missing IDs, then scrape:

```bash
.venv312/bin/python phosphosite_agent_tool.py run --protein-names-file protein_names.txt --continue-on-error
```

Resolve missing IDs, scrape, then classify the output against an expected manifest:

```bash
.venv312/bin/python phosphosite_agent_tool.py run \
  --protein-names-file protein_names.txt \
  --continue-on-error \
  --eval-manifest scraper_eval/tests/fixtures/manifest.csv \
  --eval-report scraper_eval_report.json
```

Scrape known IDs directly:

```bash
.venv312/bin/python phosphosite_agent_tool.py scrape --protein-ids-file protein_ids.txt --continue-on-error
```

The tool prints a final `SUMMARY_JSON: ...` line for agents to parse. When `--eval-manifest` is supplied, that summary includes post-scrape `evaluation` with `failure_class`, `failure_signals`, `recommended_action`, `usable_scrape_rate`, `cloudflare_likely_rate`, and `wrong_protein_rate`.

## Cloudflare-Safe Run Pattern

The scraper does not bypass Cloudflare. It is designed to reduce repeated failed traffic and preserve useful progress:

- Detect challenge pages.
- Stop aggressive retries.
- Save lookup state and partial resolved IDs after each attempt.
- Reuse cached protein IDs and completed resolved-ID scrape checkpoints.
- Back off with randomized jitter.
- Resume later from the next missing protein.

Use `--delay` as the base wait and `--delay-jitter` to randomize waits around that value:

```bash
.venv312/bin/python phosphosite_agent_tool.py run \
  --protein-names-file protein_names.txt \
  --delay 20 \
  --delay-jitter 0.5 \
  --cloudflare-cooldown 1800 \
  --max-sites-per-run 5 \
  --continue-on-error
```

With `--delay 20 --delay-jitter 0.5`, each between-request wait is randomly chosen from roughly 10 to 30 seconds. Cloudflare retry waits use randomized exponential backoff.
With `--max-sites-per-run 5`, resolved-ID scraping checkpoints and stops voluntarily after five attempted site pages. Cloudflare-blocked sites are retried after never-attempted sites.

For already-resolved protein IDs, the scraper writes:

```text
curated_protein_ids/scrape_state.json
```

That file tracks protein-page/site-discovery progress, discovered site IDs, completed site IDs, Cloudflare-blocked site IDs, failed site IDs, and output CSV paths. On the next run, completed site IDs are skipped unless you pass `--force-rescrape`.

## Validation Tool

Validate the curated ID cache and scraped CSV outputs:

```bash
.venv312/bin/python validate_scrape_outputs.py
```

The validator checks that curated IDs have expected columns/URLs, scraped CSVs are readable, output files have rows/columns, and PubMed IDs look numeric when present. It prints a final `VALIDATION_JSON: ...` line for agents to parse.

## Protein Identity Lookup

Build a lookup table that maps raw PhosphoSitePlus names from scraper outputs to canonical human gene symbols and UniProt accessions:

```bash
.venv312/bin/python protein_identity_mapper.py \
  --all-from-outputs \
  --output-root . \
  --output-csv curated_protein_ids/protein_identity_lookup.csv
```

Or resolve pasted names directly:

```bash
.venv312/bin/python protein_identity_mapper.py \
  --names Akt1 MDM2 TRAF6 DRD2 DRD3 KPNB1 \
  --output-csv curated_protein_ids/protein_identity_lookup.csv
```

The all-output mode scans every CSV under `--output-root` for `Protein`, `Downstream protein`, and `Upstream protein` columns, deduplicates every raw protein name it finds, then queries UniProt, HGNC, NCBI Gene, and Ensembl. It keeps both the raw scraped name and normalized identity fields such as `canonical_gene`, `uniprot_accession`, `recommended_name`, `aliases`, `confidence`, and `sources`.

After the lookup table exists, write enriched copies of every scraper CSV:

```bash
.venv312/bin/python protein_identity_mapper.py \
  --enrich-outputs \
  --output-root . \
  --lookup-csv curated_protein_ids/protein_identity_lookup.csv \
  --enriched-root identity_enriched_outputs
```

The enriched CSVs preserve the original columns and add canonical gene, UniProt accession, confidence, and source columns for `Protein`, `Downstream protein`, and `Upstream protein`.

## Identity-First Scraper Evaluation

Use `scraper_eval/` when you want to score actual scraper outputs against an expected manifest.
The dedicated evaluator-agent instructions live in `scraper_eval/AGENT.md`.

```bash
.venv312/bin/python -m scraper_eval.evaluate_scraper_outputs \
  --manifest scraper_eval/tests/fixtures/manifest.csv \
  --curated-ids scraper_eval/tests/fixtures/correct_ids.csv \
  --output-root scraper_eval/tests/fixtures/correct_outputs
```

This evaluator checks protein identity before PTM scoring. Wrong protein ID, organism, or UniProt/SwissProt accession is a critical identity failure and forces `final_score = 0`.
It also classifies each protein result as `VALID_SCRAPE`, `IDENTITY_FAILURE`, `CLOUDFLARE_LIKELY`, `PARTIAL_SCRAPE`, `SCHEMA_FAILURE`, or `EMPTY_OUTPUT`.

## Direct Scraper

Look up names without scraping:

```bash
python phospho_group_scraper.py --protein-names-file protein_names.txt --lookup-only --lookup-output curated_protein_ids/resolved_protein_ids.csv
```

Run a single protein ID:

```bash
python phospho_group_scraper.py --protein-id 465
```

Run a protein-ID batch:

```bash
python phospho_group_scraper.py --protein-ids-file protein_ids.txt --continue-on-error --delay 20 --delay-jitter 0.5
```

Use `protein_names.example.txt` and `protein_ids.example.txt` as templates.
