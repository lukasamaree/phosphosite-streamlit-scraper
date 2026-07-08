# AGENTS.md

## Scraper Workflow

- This project scrapes PhosphoSitePlus proteins with `phospho_group_scraper.py`.
- Codex should reason about repo state before running a tool. Do not hide this choice behind a deterministic planner script.
- First inspect the curated ID cache (`curated_protein_ids/resolved_protein_ids.csv`) and lookup checkpoint (`curated_protein_ids/lookup_state.json`) when they exist.
- Then inspect existing protein output folders/CSVs for requested proteins before deciding to scrape.
- Choose the smallest safe command:
  - Resolve IDs only when protein names lack saved human PhosphoSitePlus IDs.
  - Scrape known IDs directly when curator IDs already exist.
  - Run the full workflow only when both ID resolution and scraping are needed.
  - Rerun only failed or missing proteins unless the user asks for a full refresh.
- Prefer the reusable tool wrapper `phosphosite_agent_tool.py` for execution after Codex has made the decision.
- After scraping, validate outputs with `.venv312/bin/python validate_scrape_outputs.py` and parse the final `VALIDATION_JSON: ...` line.
- For deterministic scraper-quality scoring, use the separate evaluator-agent instructions in `scraper_eval/AGENT.md`. It gates PTM scoring behind protein identity and prints `EVAL_JSON: ...`.
- Resolve and save IDs for a worklist with `.venv312/bin/python phosphosite_agent_tool.py resolve --protein-names-file protein_names.txt`.
- Scrape from saved IDs, resolving missing names first, with `.venv312/bin/python phosphosite_agent_tool.py run --protein-names-file protein_names.txt --continue-on-error`.
- Scrape already-known IDs with `.venv312/bin/python phosphosite_agent_tool.py scrape --protein-ids-file protein_ids.txt --continue-on-error`.
- After running a tool, parse the final `SUMMARY_JSON: ...` line, inspect generated CSVs or failures, and adapt the next command from that evidence.
- After validation, treat `status=failed` as a blocker. Treat `passed_with_warnings` as usable but worth reviewing.
- Treat any `critical_identity_failure` from `scraper_eval` as a blocker, even if PTM rows are present.
- The reusable tool and dashboard both write resolved IDs to `curated_protein_ids/`; generated scrape output folders and site CSVs should remain uncommitted.
- Curators may maintain a local `protein_names.txt` file with one protein name/symbol per line. This is the preferred workflow.
- The scraper looks up protein names by searching PhosphoSitePlus, choosing a human result row, navigating to the human protein page, and extracting the ID from the final resolved `proteinAction.action?id=...` URL. Treat the final URL as the source of truth.
- Curator IDs, when used directly, are the `id=` values from protein URLs such as `https://www.phosphosite.org/proteinAction.action?id=582&showAllSites=true`.
- Curators may also maintain a local `protein_ids.txt` file with one protein ID per line. Do not commit real curator worklists unless explicitly requested.
- Look up names without scraping with `python phospho_group_scraper.py --protein-names-file protein_names.txt --lookup-only --lookup-output curated_protein_ids/resolved_protein_ids.csv`.
- For name lookup, use only human protein results and use the resolved `proteinAction.action?id=...` URL as the authority. Do not infer IDs from visible text alone.
- Run a protein-name worklist with `python phospho_group_scraper.py --protein-names-file protein_names.txt --continue-on-error --delay 2`.
- Run a single protein ID with `python phospho_group_scraper.py --protein-id 582`.
- Run a curator worklist with `python phospho_group_scraper.py --protein-ids-file protein_ids.txt --continue-on-error --delay 2`.
- Direct `siteAction` IDs are supported only for advanced troubleshooting with `--site-id`, `--site-ids`, or `--site-ids-file`.
- Keep request rates conservative. Do not reduce `--delay` or add concurrency unless explicitly requested.
- Generated protein folders and CSV outputs are scrape artifacts. Do not commit them unless the user asks.
- If a run fails, classify the failure as dependency, browser installation, network/access, invalid protein ID, no site links found, or scraper selector/data-shape breakage before changing code.
- After changing scraper behavior, verify the CLI with `python phospho_group_scraper.py --help` and, when network/browser access is available, a one-protein scrape.
