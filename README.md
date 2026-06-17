# PhosphoSitePlus Streamlit Scraper

Streamlit dashboard and reusable command-line tools for resolving human PhosphoSitePlus protein URL IDs and scraping phosphorylation/site data.

## What It Does

- Accepts protein symbols such as `AKT1`, `TP53`, `EGFR`, and `MAPK1`.
- Searches PhosphoSitePlus and selects the human protein page.
- Reads the protein ID from the final `proteinAction.action?id=...` URL.
- Saves curated protein IDs in `curated_protein_ids/`.
- Reuses saved IDs so known proteins do not need to be looked up again.
- Scrapes resolved protein IDs with conservative delays and live logs.

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

Scrape known IDs directly:

```bash
.venv312/bin/python phosphosite_agent_tool.py scrape --protein-ids-file protein_ids.txt --continue-on-error
```

The tool prints a final `SUMMARY_JSON: ...` line for agents to parse.

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
python phospho_group_scraper.py --protein-ids-file protein_ids.txt --continue-on-error --delay 2
```

Use `protein_names.example.txt` and `protein_ids.example.txt` as templates.
