import io
import json
import os
import shlex
import subprocess
import sys
import tempfile
import time
import zipfile
from pathlib import Path

import pandas as pd
import streamlit as st


ROOT = Path(__file__).resolve().parent
SCRAPER = ROOT / "phospho_group_scraper.py"
AGENTIC_WORKFLOW = ROOT / "agentic_phospho_workflow.py"
IDENTITY_MAPPER = ROOT / "protein_identity_mapper.py"
CURATION_DIR = ROOT / "curated_protein_ids"
LOOKUP_CSV = CURATION_DIR / "resolved_protein_ids.csv"
IDENTITY_LOOKUP_CSV = CURATION_DIR / "protein_identity_lookup.csv"
IDS_TXT = CURATION_DIR / "protein_ids.txt"
AGENTIC_STATE_JSON = CURATION_DIR / "lookup_state.json"
SCRAPE_STATE_JSON = CURATION_DIR / "scrape_state.json"

CURATION_DIR.mkdir(exist_ok=True)


st.set_page_config(
    page_title="PhosphoSitePlus Scraper",
    layout="wide",
    initial_sidebar_state="expanded",
)


@st.cache_resource(show_spinner="Checking Playwright browser...")
def ensure_playwright_browser():
    if os.environ.get("PHOSPHOSITE_SKIP_PLAYWRIGHT_CHECK") == "1":
        return "skipped"

    install_command = [sys.executable, "-m", "playwright", "install", "chromium"]
    completed = subprocess.run(
        install_command,
        cwd=ROOT,
        text=True,
        capture_output=True,
        timeout=300,
    )
    if completed.returncode != 0:
        output = "\n".join(part for part in [completed.stdout, completed.stderr] if part)
        raise RuntimeError(f"Playwright Chromium install/check failed:\n{output}")
    return "installed/verified"


def parse_protein_text(text):
    names = []
    normalized = text.replace("[", "\n").replace("]", "\n").replace(",", "\n")
    for line in normalized.splitlines():
        value = line.split("#", 1)[0].strip()
        value = value.strip("'\" \t")
        if value:
            names.append(value)
    return names


def unique_names(names):
    seen = set()
    unique = []
    for name in names:
        key = name.strip().upper()
        if key and key not in seen:
            seen.add(key)
            unique.append(name.strip())
    return unique


def run_scraper_command(args, timeout=None):
    command = [sys.executable, "-u", "-B", str(SCRAPER), *args]
    completed = subprocess.run(
        command,
        cwd=ROOT,
        text=True,
        capture_output=True,
        timeout=timeout,
    )
    output = "\n".join(part for part in [completed.stdout, completed.stderr] if part)
    return completed.returncode, output


def run_agentic_command(args, timeout=None):
    command = [sys.executable, "-u", "-B", str(AGENTIC_WORKFLOW), *args]
    completed = subprocess.run(
        command,
        cwd=ROOT,
        text=True,
        capture_output=True,
        timeout=timeout,
    )
    output = "\n".join(part for part in [completed.stdout, completed.stderr] if part)
    return completed.returncode, output


def stream_command(script_path, args, label, summary_lines=None):
    command = [sys.executable, "-u", "-B", str(script_path), *args]
    status_box = st.empty()
    log_box = st.empty()
    lines = []
    start_time = time.monotonic()

    def add_log(message):
        elapsed = time.monotonic() - start_time
        line = f"[{elapsed:7.1f}s] {message}"
        lines.append(line)
        status_box.info(message)
        log_box.code("\n".join(lines[-120:]), language="text")

    add_log(f"{label}: preparing subprocess")
    for summary_line in summary_lines or []:
        add_log(summary_line)
    add_log("Checking Playwright Chromium before scraper run")
    playwright_status = ensure_playwright_browser()
    add_log(f"Playwright Chromium: {playwright_status}")
    add_log("Command: " + " ".join(shlex.quote(part) for part in command))
    process = subprocess.Popen(
        command,
        cwd=ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        bufsize=1,
    )

    assert process.stdout is not None
    for line in process.stdout:
        clean_line = line.rstrip()
        if not clean_line:
            continue
        add_log(clean_line)

    return_code = process.wait()
    output = "\n".join(lines)
    if return_code == 0:
        status_box.success(f"{label}: finished in {time.monotonic() - start_time:.1f}s")
    else:
        status_box.warning(f"{label}: exited with code {return_code} after {time.monotonic() - start_time:.1f}s")
    return return_code, output


def read_lookup_csv(path=LOOKUP_CSV):
    if not path.exists():
        return pd.DataFrame(columns=["protein_name", "protein_id", "url", "organism", "source"])
    return pd.read_csv(path)


def read_agentic_state(path=AGENTIC_STATE_JSON):
    if not path.exists():
        return {"runs": {}, "updated_at": None}
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def resolved_rows_for_proteins(proteins, lookup_df):
    if lookup_df.empty:
        return lookup_df, proteins
    if not proteins:
        return lookup_df, []

    requested_keys = {protein.strip().upper() for protein in proteins}
    df = lookup_df.copy()
    df["_protein_key"] = df["protein_name"].astype(str).str.strip().str.upper()
    found = df[df["_protein_key"].isin(requested_keys)].drop(columns=["_protein_key"])
    found_keys = set(found["protein_name"].astype(str).str.strip().str.upper())
    missing = [protein for protein in proteins if protein.strip().upper() not in found_keys]
    return found, missing


def write_ids_file(ids, path=IDS_TXT):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        for protein_id in ids:
            handle.write(f"{int(protein_id)}\n")


@st.cache_data(ttl=30)
def output_files():
    ignored = {".git", ".venv", "__pycache__"}
    files = []
    for path in ROOT.rglob("*.csv"):
        relative = path.relative_to(ROOT)
        if any(part in ignored for part in relative.parts):
            continue
        if str(relative).startswith("Home >"):
            continue
        files.append(path)
    return sorted(files, key=lambda item: str(item.relative_to(ROOT)).lower())


@st.cache_data(ttl=30)
def summarize_outputs(files):
    rows = []
    for path in files:
        try:
            with open(path, "r", encoding="utf-8", errors="replace") as handle:
                row_count = max(sum(1 for _ in handle) - 1, 0)
            header = pd.read_csv(path, nrows=0)
            rows.append(
                {
                    "protein": path.parent.name,
                    "site": path.stem.split("_", 1)[0],
                    "file": str(path.relative_to(ROOT)),
                    "rows": row_count,
                    "columns": len(header.columns),
                    "size_kb": round(path.stat().st_size / 1024, 1),
                }
            )
        except Exception as exc:
            rows.append(
                {
                    "protein": path.parent.name,
                    "site": path.stem,
                    "file": str(path.relative_to(ROOT)),
                    "rows": None,
                    "columns": None,
                    "size_kb": round(path.stat().st_size / 1024, 1),
                    "error": str(exc),
                }
            )
    return pd.DataFrame(rows)


def dataframe_download(df, filename, label):
    buffer = io.StringIO()
    df.to_csv(buffer, index=False)
    st.download_button(label, buffer.getvalue(), file_name=filename, mime="text/csv")


def zip_files(files):
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for path in files:
            archive.write(path, arcname=str(path.relative_to(ROOT)))
    buffer.seek(0)
    return buffer.getvalue()


st.title("PhosphoSitePlus Scraper Dashboard")
st.caption(
    "Human-only workflow: search protein name, navigate to the human protein page, "
    "then read the ID from the final proteinAction URL."
)

playwright_status = "checked before scraper runs"

if "protein_text" not in st.session_state:
    st.session_state.protein_text = "AKT1\nTP53\nEGFR\nMAPK1"

with st.sidebar:
    st.subheader("Workflow")
    delay = st.number_input("Delay between requests", min_value=0.0, max_value=30.0, value=5.0, step=0.5)
    delay_jitter = st.slider("Delay jitter", min_value=0.0, max_value=1.0, value=0.5, step=0.05)
    gentle_mode = st.toggle("Gentle mode", value=True)
    max_sites_per_run = st.number_input("Max sites per run", min_value=0, max_value=100, value=5, step=1)
    attempts = st.number_input("Agent retry attempts", min_value=1, max_value=10, value=3, step=1)
    cloudflare_cooldown = st.number_input(
        "Cloudflare cooldown",
        min_value=10.0,
        max_value=1800.0,
        value=900.0,
        step=30.0,
    )
    continue_on_error = st.toggle("Continue on errors", value=True)
    reuse_saved_ids = st.toggle("Reuse saved protein IDs", value=True)
    resolve_missing_before_scrape = st.toggle("Resolve missing IDs before scrape", value=True)
    st.divider()
    st.caption("Outputs are written under this repository.")
    st.code(str(ROOT), language="text")
    st.caption(f"Playwright browser: {playwright_status}")

effective_delay = max(delay, 20.0) if gentle_mode else delay
effective_cloudflare_cooldown = max(cloudflare_cooldown, 1800.0) if gentle_mode else cloudflare_cooldown
effective_max_sites_per_run = max_sites_per_run if gentle_mode else max_sites_per_run

run_tab, lookup_tab, identity_tab, outputs_tab = st.tabs(["Run", "Resolved IDs", "Identity Lookup", "Outputs"])

with run_tab:
    left, right = st.columns([0.9, 1.1], gap="large")

    with left:
        st.subheader("Protein Worklist")
        protein_text = st.text_area(
            "Protein symbols",
            key="protein_text",
            height=180,
            placeholder="AKT1\nTP53\nEGFR\nMAPK1",
        )
        proteins = unique_names(parse_protein_text(protein_text))
        st.caption(f"{len(proteins)} protein symbol(s) ready")

        saved_lookup_df = read_lookup_csv()
        saved_for_input, missing_for_input = resolved_rows_for_proteins(proteins, saved_lookup_df)
        if proteins and reuse_saved_ids:
            st.caption(f"{len(saved_for_input)} saved ID(s), {len(missing_for_input)} missing for this worklist")

        lookup_clicked = st.button("Resolve Human URL IDs", type="primary", use_container_width=True)
        agentic_lookup_clicked = st.button("Agentic Resolve IDs", use_container_width=True)
        scrape_clicked = st.button("Scrape Resolved IDs", use_container_width=True)

    with right:
        st.subheader("Run Log")
        if lookup_clicked:
            if not proteins:
                st.error("Add at least one protein symbol.")
            else:
                with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False, encoding="utf-8") as handle:
                    for protein in proteins:
                        handle.write(f"{protein}\n")
                    names_file = handle.name

                args = [
                    "--protein-names-file",
                    names_file,
                    "--lookup-only",
                    "--lookup-output",
                    str(LOOKUP_CSV),
                    "--delay",
                    str(effective_delay),
                    "--delay-jitter",
                    str(delay_jitter),
                ]
                if continue_on_error:
                    args.append("--continue-on-error")

                summary_lines = [
                    f"Input: {len(proteins)} protein symbol(s)",
                    f"Output CSV: {LOOKUP_CSV}",
                    f"Delay between proteins: {effective_delay}s +/- {delay_jitter:.0%}",
                    "Workflow: search name -> choose human result -> read id= from final URL",
                ]
                with st.spinner("Searching names, opening human protein pages, and reading final URL IDs..."):
                    code, output = stream_command(SCRAPER, args, "Human URL lookup", summary_lines)

                if code == 0:
                    df = read_lookup_csv()
                    st.success(f"Resolved {len(df)} protein ID(s).")
                    if not df.empty:
                        write_ids_file(df["protein_id"].dropna().astype(int).tolist())
                        st.dataframe(df, use_container_width=True, hide_index=True)
                else:
                    st.error("Lookup command failed.")

        if agentic_lookup_clicked:
            if not proteins:
                st.error("Add at least one protein symbol.")
            else:
                with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False, encoding="utf-8") as handle:
                    for protein in proteins:
                        handle.write(f"{protein}\n")
                    names_file = handle.name

                args = [
                    "--protein-names-file",
                    names_file,
                    "--output-csv",
                    str(LOOKUP_CSV),
                    "--state-json",
                    str(AGENTIC_STATE_JSON),
                    "--attempts",
                    str(attempts),
                    "--delay",
                    str(effective_delay),
                    "--delay-jitter",
                    str(delay_jitter),
                    "--cloudflare-cooldown",
                    str(effective_cloudflare_cooldown),
                ]

                summary_lines = [
                    f"Input: {len(proteins)} protein symbol(s)",
                    f"Attempts per protein: {attempts}",
                    f"Delay/backoff base: {effective_delay}s +/- {delay_jitter:.0%}",
                    f"Cloudflare cooldown: {effective_cloudflare_cooldown}s",
                    f"Checkpoint JSON: {AGENTIC_STATE_JSON}",
                    f"Partial CSV: {LOOKUP_CSV}",
                    "Workflow: resume resolved proteins, retry failures, checkpoint after every attempt",
                ]
                with st.spinner("Agent is resolving IDs with checkpointing, retries, and backoff..."):
                    code, output = stream_command(AGENTIC_WORKFLOW, args, "Agentic lookup", summary_lines)

                df = read_lookup_csv()
                if not df.empty:
                    write_ids_file(df["protein_id"].dropna().astype(int).tolist())
                    st.dataframe(df, use_container_width=True, hide_index=True)
                if code == 0:
                    st.success(f"Agentic lookup finished with {len(df)} resolved protein ID(s).")
                else:
                    st.warning("Agentic lookup finished with errors. Partial results were checkpointed.")

        if scrape_clicked:
            lookup_df = read_lookup_csv()
            requested_proteins = proteins
            usable_lookup_df, missing_proteins = resolved_rows_for_proteins(requested_proteins, lookup_df)

            if requested_proteins and resolve_missing_before_scrape and missing_proteins:
                st.info(f"Resolving {len(missing_proteins)} missing protein ID(s) before scraping.")
                with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False, encoding="utf-8") as handle:
                    for protein in missing_proteins:
                        handle.write(f"{protein}\n")
                    names_file = handle.name

                lookup_args = [
                    "--protein-names-file",
                    names_file,
                    "--output-csv",
                    str(LOOKUP_CSV),
                    "--state-json",
                    str(AGENTIC_STATE_JSON),
                    "--attempts",
                    str(attempts),
                    "--delay",
                    str(effective_delay),
                    "--delay-jitter",
                    str(delay_jitter),
                    "--cloudflare-cooldown",
                    str(effective_cloudflare_cooldown),
                ]
                summary_lines = [
                    f"Saved IDs available: {len(usable_lookup_df)}",
                    f"Missing IDs to resolve now: {len(missing_proteins)}",
                    f"Missing proteins: {', '.join(missing_proteins)}",
                    f"Delay/backoff base: {effective_delay}s +/- {delay_jitter:.0%}",
                    f"Cloudflare cooldown: {effective_cloudflare_cooldown}s",
                    f"Checkpoint JSON: {AGENTIC_STATE_JSON}",
                    "Workflow: use saved IDs first, resolve only missing proteins, then scrape",
                ]
                with st.spinner("Resolving missing protein IDs before scraping..."):
                    lookup_code, lookup_output = stream_command(
                        AGENTIC_WORKFLOW,
                        lookup_args,
                        "Resolve missing IDs",
                        summary_lines,
                    )

                lookup_df = read_lookup_csv()
                usable_lookup_df, missing_proteins = resolved_rows_for_proteins(requested_proteins, lookup_df)
                if lookup_code != 0:
                    st.warning("Missing-ID lookup finished with errors. Scraping any IDs that were resolved.")

            elif requested_proteins and missing_proteins and not resolve_missing_before_scrape:
                st.warning(f"Skipping {len(missing_proteins)} protein(s) without saved IDs: {', '.join(missing_proteins)}")

            if requested_proteins:
                lookup_df = usable_lookup_df

            if lookup_df.empty:
                st.error("No resolved IDs available. Resolve IDs first, or enable missing-ID lookup before scraping.")
            else:
                ids = lookup_df["protein_id"].dropna().astype(int).tolist()
                write_ids_file(ids)
                args = [
                    "--protein-ids-file",
                    str(IDS_TXT),
                    "--delay",
                    str(effective_delay),
                    "--delay-jitter",
                    str(delay_jitter),
                    "--cloudflare-cooldown",
                    str(effective_cloudflare_cooldown),
                    "--max-sites-per-run",
                    str(effective_max_sites_per_run),
                    "--scrape-state",
                    str(SCRAPE_STATE_JSON),
                ]
                if continue_on_error:
                    args.append("--continue-on-error")

                summary_lines = [
                    f"Input: {len(ids)} resolved protein ID(s)",
                    f"ID file: {IDS_TXT}",
                    f"Saved ID CSV: {LOOKUP_CSV}",
                    f"Scrape checkpoint JSON: {SCRAPE_STATE_JSON}",
                    f"Gentle mode: {'on' if gentle_mode else 'off'}",
                    f"Max sites per run: {effective_max_sites_per_run or 'unlimited'}",
                    f"Cloudflare circuit-breaker cooldown: {effective_cloudflare_cooldown}s",
                    f"Missing proteins skipped: {len(missing_proteins) if requested_proteins else 0}",
                    f"Delay between IDs/sites: {effective_delay}s +/- {delay_jitter:.0%}",
                    "Workflow: open protein page -> collect human siteAction links -> scrape each site -> write CSVs",
                ]
                with st.spinner(f"Scraping {len(ids)} protein ID(s)..."):
                    code, output = stream_command(SCRAPER, args, "Protein scrape", summary_lines)

                if code == 0:
                    st.success("Scrape completed.")
                else:
                    st.warning("Scrape finished with errors. Check the log for failed IDs.")

with lookup_tab:
    st.subheader("Resolved Protein IDs")
    lookup_df = read_lookup_csv()

    if lookup_df.empty:
        st.info("Run human URL ID lookup from the Run tab to populate this table.")
    else:
        st.dataframe(lookup_df, use_container_width=True, hide_index=True)
        dataframe_download(lookup_df, "resolved_protein_ids.csv", "Download ID CSV")

    st.divider()
    st.subheader("Manual IDs")
    manual_ids = st.text_area(
        "Paste PhosphoSitePlus protein IDs",
        value="\n".join(str(value) for value in lookup_df.get("protein_id", pd.Series(dtype=int)).dropna().astype(int).tolist()),
        height=140,
        placeholder="570\n465\n592\n832",
    )
    if st.button("Save Manual ID File"):
        ids = []
        for line in manual_ids.splitlines():
            value = line.strip()
            if value:
                ids.append(int(value))
        write_ids_file(ids)
        st.success(f"Saved {len(ids)} ID(s) to {IDS_TXT.name}.")

with identity_tab:
    st.subheader("Protein Identity Lookup")
    st.caption("Map raw PhosphoSite protein names to canonical human gene symbols and UniProt accessions.")

    identity_names_text = st.text_area(
        "Optional raw protein names",
        height=120,
        placeholder="Akt1\nMDM2\nTRAF6\nDRD2\nDRD3\nKPNB1",
    )
    col_a, col_b = st.columns(2)
    build_from_outputs = col_a.button("Build From Output CSVs", use_container_width=True)
    build_from_names = col_b.button("Build From Pasted Names", use_container_width=True)

    if build_from_outputs or build_from_names:
        args = [
            "--output-csv",
            str(IDENTITY_LOOKUP_CSV),
            "--cache-json",
            str(CURATION_DIR / "protein_identity_lookup_cache.json"),
        ]
        summary_lines = [
            f"Output CSV: {IDENTITY_LOOKUP_CSV}",
            "Sources: UniProt, HGNC, NCBI Gene, Ensembl",
            "Columns scanned: Protein, Downstream protein, Upstream protein",
        ]
        if build_from_names:
            names = unique_names(parse_protein_text(identity_names_text))
            if not names:
                st.error("Paste at least one raw protein name.")
                args = None
            else:
                args.extend(["--names", *names])
                summary_lines.insert(0, f"Input names: {len(names)}")
        if args is not None:
            with st.spinner("Building canonical gene/UniProt lookup table..."):
                code, output = stream_command(IDENTITY_MAPPER, args, "Protein identity lookup", summary_lines)
            if code == 0:
                st.success("Protein identity lookup table built.")
            else:
                st.warning("Protein identity lookup finished with errors. Check the log.")

    if IDENTITY_LOOKUP_CSV.exists():
        identity_df = pd.read_csv(IDENTITY_LOOKUP_CSV)
        st.dataframe(identity_df, use_container_width=True, hide_index=True)
        dataframe_download(identity_df, "protein_identity_lookup.csv", "Download Identity Lookup CSV")
    else:
        st.info("Build the identity lookup table from outputs or pasted raw names.")

with outputs_tab:
    st.subheader("Scraped Outputs")
    if st.button("Scan Output Folders", use_container_width=True):
        st.session_state.outputs_scanned = True

    if not st.session_state.get("outputs_scanned"):
        st.info("Click Scan Output Folders to load CSV summaries.")
    else:
        files = output_files()
        summary = summarize_outputs(files)

        metric_cols = st.columns(4)
        metric_cols[0].metric("CSV files", len(files))
        metric_cols[1].metric("Proteins", summary["protein"].nunique() if not summary.empty else 0)
        metric_cols[2].metric("Total rows", int(summary["rows"].fillna(0).sum()) if "rows" in summary else 0)
        metric_cols[3].metric("Total size", f"{summary['size_kb'].sum():.1f} KB" if not summary.empty else "0 KB")

        if summary.empty:
            st.info("No output CSV files found yet.")
        else:
            protein_options = ["All", *sorted(summary["protein"].dropna().unique())]
            selected_protein = st.selectbox("Protein", protein_options)
            filtered = summary if selected_protein == "All" else summary[summary["protein"] == selected_protein]
            filtered_files = [ROOT / file_path for file_path in filtered["file"].tolist()]

            st.dataframe(filtered, use_container_width=True, hide_index=True)
            dataframe_download(filtered, "scrape_output_summary.csv", "Download Output Summary")
            st.download_button(
                "Download Output ZIP",
                zip_files(filtered_files),
                file_name="phosphosite_scrape_outputs.zip",
                mime="application/zip",
            )

            counts = (
                filtered.groupby("protein", dropna=False)["file"]
                .count()
                .sort_values(ascending=False)
                .reset_index(name="csv_files")
            )
            st.dataframe(counts, use_container_width=True, hide_index=True)

            selected_file = st.selectbox("Preview CSV", filtered["file"].tolist())
            preview_path = ROOT / selected_file
            preview_df = pd.read_csv(preview_path)
            st.dataframe(preview_df, use_container_width=True)
            dataframe_download(preview_df, preview_path.name, "Download Selected CSV")
