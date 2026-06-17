# Deploying The Streamlit App

Use Streamlit Community Cloud with:

- Repository: `lukasamaree/phosphosite-streamlit-scraper`
- Branch: `main`
- Main file path: `streamlit_phospho_scraper.py`
- Python runtime: `python-3.12`

Generated scrape outputs are intentionally ignored. The lightweight curated ID cache is kept in `curated_protein_ids/` so resolved protein IDs can be committed when you want to preserve them.

For local development:

```bash
python3.12 -m venv .venv312
.venv312/bin/python -m pip install -r requirements.txt
.venv312/bin/python -m streamlit run streamlit_phospho_scraper.py
```
