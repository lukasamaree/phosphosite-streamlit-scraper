import asyncio
import argparse
import csv
import re
import os
import random
from pathlib import Path
from urllib.parse import quote_plus


DEFAULT_CLOUDFLARE_RETRIES = 3
DEFAULT_CLOUDFLARE_WAIT_SECONDS = 10
DEFAULT_DELAY_JITTER = 0.5
ROOT = Path(__file__).resolve().parent
BROWSER_STATE_PATH = ROOT / ".phosphosite_browser_state.json"


def cloudflare_wait_seconds(base_wait, attempt):
    exponential_wait = base_wait * (2 ** (attempt - 1))
    return random.uniform(exponential_wait, exponential_wait * 1.75)


def jittered_delay_seconds(base_delay, jitter=DEFAULT_DELAY_JITTER):
    if base_delay <= 0:
        return 0
    jitter = max(0.0, jitter)
    lower = max(0.0, base_delay * (1 - jitter))
    upper = base_delay * (1 + jitter)
    return random.uniform(lower, upper)


async def sleep_between_requests(base_delay, jitter, label):
    wait_for = jittered_delay_seconds(base_delay, jitter)
    if wait_for <= 0:
        return
    print(
        f"WAIT: sleeping {wait_for:.1f}s before next {label} "
        f"(base={base_delay:.1f}s, jitter={jitter:.2f})",
        flush=True,
    )
    await asyncio.sleep(wait_for)


async def is_cloudflare_challenge_page(page):
    try:
        title = (await page.title()).strip().lower()
    except Exception:
        title = ""
    if "just a moment" in title or "checking your browser" in title:
        return True

    try:
        body_text = (await page.locator("body").inner_text(timeout=2000)).lower()
    except Exception:
        body_text = ""
    challenge_markers = [
        "cloudflare",
        "verify you are human",
        "checking if the site connection is secure",
        "enable javascript and cookies",
    ]
    return any(marker in body_text for marker in challenge_markers)


async def goto_with_cloudflare_retry(
    page,
    url,
    label=None,
    retries=DEFAULT_CLOUDFLARE_RETRIES,
    wait_seconds=DEFAULT_CLOUDFLARE_WAIT_SECONDS,
):
    label = label or url
    last_title = ""
    for attempt in range(1, retries + 1):
        await page.goto(url, wait_until="domcontentloaded", timeout=60000)
        try:
            last_title = await page.title()
        except Exception:
            last_title = ""

        if not await is_cloudflare_challenge_page(page):
            return

        print(
            f"CLOUDFLARE: challenge detected while loading {label}; "
            f"attempt {attempt}/{retries}; waiting {wait_seconds}s",
            flush=True,
        )
        wait_for = cloudflare_wait_seconds(wait_seconds, attempt)
        print(f"CLOUDFLARE: waiting {wait_for:.1f}s before checking {label} again", flush=True)
        await page.wait_for_timeout(wait_for * 1000)
        if not await is_cloudflare_challenge_page(page):
            print(f"CLOUDFLARE: challenge cleared for {label}", flush=True)
            return

    raise RuntimeError(
        f"Cloudflare challenge persisted while loading {label} "
        f"after {retries} attempt(s); last title='{last_title}'."
    )


async def new_browser_context(playwright):
    browser = await playwright.chromium.launch(headless=True)
    context_args = {
        "viewport": {"width": 1365, "height": 900},
        "user_agent": (
            "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
            "AppleWebKit/537.36 (KHTML, like Gecko) "
            "Chrome/125.0.0.0 Safari/537.36"
        ),
    }
    if BROWSER_STATE_PATH.exists() and BROWSER_STATE_PATH.stat().st_size > 0:
        context_args["storage_state"] = str(BROWSER_STATE_PATH)
        print(f"SESSION: reusing browser storage from {BROWSER_STATE_PATH}", flush=True)
    context = await browser.new_context(**context_args)
    page = await context.new_page()
    return browser, context, page


async def close_browser_context(browser, context):
    try:
        await context.storage_state(path=str(BROWSER_STATE_PATH))
        print(f"SESSION: saved browser storage to {BROWSER_STATE_PATH}", flush=True)
    except Exception as exc:
        print(f"WARNING: could not save browser storage: {exc}", flush=True)
    await browser.close()



# Protein IDs come from PhosphoSitePlus protein URLs:
# https://www.phosphosite.org/proteinAction.action?id=582&showAllSites=true

def first_webscraper(page):
    """
    Extracts the protein name and phosphosite from the header, removing (human) from the protein name.
    Returns (amino_acid, protein_name)
    """
    import re
    async def inner():
        amino_acid = None
        protein_name = None
        header_div = await page.query_selector('#titleMainHeader')
        if header_div:
            header_text = await header_div.inner_text()
            # Example: "Phosphorylation Site Page: > Thr160 - CDK2 (human)"
            match = re.search(r'>\s*([A-Za-z0-9]+)\s*-\s*([A-Za-z0-9_\-]+)\s*\(human\)', header_text)
            if match:
                amino_acid = match.group(1)
                protein_name = match.group(2)
            else:
                amino_acid = header_text
                protein_name = header_text
        if protein_name:
            protein_name = re.sub(r'\(human\)', '', protein_name, flags=re.IGNORECASE).replace('(', '').replace(')', '').strip()
        return amino_acid, protein_name
    return inner


def safe_filename_component(value):
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value or "unknown").strip())
    return cleaned.strip("_") or "unknown"

async def upstream_scraper(page):
    """
    Scrapes the Upstream Regulation table and returns a dictionary with keys:
    'Regulatory protein', 'Putative in vivo kinases', 'Kinases in vitro', 'Phosphatases in vitro'.
    Each value is a string containing the scraped text for that field.
    """
    result = {
        'Regulatory protein': '',
        'Putative in vivo kinases': '',
        'Kinases in vitro': '',
        'Phosphatases in vitro': ''
    }
    # Find the Upstream Regulation table by its <th> text
    tables = await page.query_selector_all('table')
    for table in tables:
        ths = await table.query_selector_all('th')
        for th in ths:
            th_text = (await th.inner_text()).strip().lower()
            if 'upstream regulation' in th_text:
                # This is the correct table
                trs = await table.query_selector_all('tr')
                for tr in trs:
                    tds = await tr.query_selector_all('td')
                    if len(tds) >= 2:
                        field = (await tds[0].inner_text()).strip().lower()
                        value = (await tds[1].inner_text()).strip()
                        if 'regulatory protein' in field:
                            result['Regulatory protein'] = value
                        elif 'putative in vivo kinases' in field:
                            result['Putative in vivo kinases'] = value
                        elif 'kinases, in vitro' in field:
                            result['Kinases in vitro'] = value
                        elif 'phosphatases, in vitro' in field:
                            result['Phosphatases in vitro'] = value
                return result
    return result

async def downstream_scraper(page, protein_name):
    """
    Scrapes the Downstream Regulation table and returns a dictionary with keys:
    'Effects of modification on {protein_name}',
    'Effects of modification on biological processes',
    'Induce interaction with:',
    'Inhibit interaction with:'.
    Each value is a string containing the scraped text for that field.
    """
    result = {
        f'Effects of modification on {protein_name}': '',
        'Effects of modification on biological processes': '',
        'Induce interaction with:': '',
        'Inhibit interaction with:': ''
    }
    # Find the Downstream Regulation table by its <th> text
    tables = await page.query_selector_all('table')
    for table in tables:
        ths = await table.query_selector_all('th')
        for th in ths:
            th_text = (await th.inner_text()).strip().lower()
            if 'downstream regulation' in th_text:
                # This is the correct table
                trs = await table.query_selector_all('tr')
                for tr in trs:
                    tds = await tr.query_selector_all('td')
                    if len(tds) >= 2:
                        field = (await tds[0].inner_text()).strip().lower()
                        value = (await tds[1].inner_text()).strip()
                        field_clean = field.strip().rstrip(':').strip()
                        if field_clean == f'effects of modification on {protein_name.lower()}':
                            result[f'Effects of modification on {protein_name}'] = value
                        elif field_clean == 'effects of modification on biological processes':
                            result['Effects of modification on biological processes'] = value
                        elif 'induce interaction with' in field:
                            result['Induce interaction with:'] = value
                        elif 'inhibit interaction with' in field:
                            result['Inhibit interaction with:'] = value
                return result
    return result

async def references_scraper(page):
    """
    Scrapes the References table and returns a list of dicts with keys 'Reference Number' and 'PubMed ID'.
    """
    results = []
    # Find the References table by its <th> text
    tables = await page.query_selector_all('table')
    for table in tables:
        ths = await table.query_selector_all('th')
        for th in ths:
            th_text = (await th.inner_text()).strip().lower()
            if th_text == 'references':
                # This is the correct table
                trs = await table.query_selector_all('tr')
                for tr in trs:
                    tds = await tr.query_selector_all('td')
                    if len(tds) >= 2:
                        # Reference number is in the first column
                        ref_number = (await tds[0].inner_text()).strip()
                        # PubMed ID is the first number in the second column (bold red)
                        pubmed_match = re.search(r'(\d{7,})', await tds[1].inner_text())
                        pubmed_id = pubmed_match.group(1) if pubmed_match else None
                        results.append({
                            'Reference Number': ref_number,
                            'PubMed ID': pubmed_id
                        })
                return results
    return results

async def main(site_id):
    from playwright.async_api import async_playwright

    url = f"https://www.phosphosite.org/siteAction.action?id={site_id}"
    async with async_playwright() as p:
        browser, context, page = await new_browser_context(p)
        await goto_with_cloudflare_retry(page, url, f"siteAction ID {site_id}")
        amino_acid, protein_name = await first_webscraper(page)()
        # Save upstream scraper result, exploded with entity, organism, references
        upstream_data = await upstream_scraper(page)
        print(f"DEBUG: Upstream data scraped: {upstream_data}")
        upstream_rows = []
        entity_pattern = re.compile(r'([A-Za-z0-9\-_,\[\] ]+?)\s*\((human|mouse)\)\s*\(([^\)]+)\)')
        for field, value in upstream_data.items():
            print(f"DEBUG: Processing field '{field}' with value: '{value}'")
            # Find all matches like NAME (organism) ( numbers )
            for match in entity_pattern.finditer(value):
                entity = match.group(1).strip()
                organism = match.group(2).strip()
                references = [int(ref.strip()) for ref in match.group(3).split(',') if ref.strip().isdigit()]
                upstream_rows.append({
                    'Upstream regulation': field,
                    'Upstream protein': entity.lstrip(', ').strip(),
                    'Organism': organism,
                    'References': references,
                    'Amino Acid': amino_acid,
                    'Protein': protein_name
                })
                print(f"DEBUG: Added row - Entity: {entity}, Organism: {organism}, References: {references}")
        print(f"DEBUG: Upstream row count: {len(upstream_rows)}")

        # Save downstream scraper result, fully exploded
        downstream_data = await downstream_scraper(page, protein_name)
        downstream_rows = []
        # Regex for effect phrase and references
        effect_pattern = re.compile(r'([^\(]+?)\s*\(([^\)]+)\)')
        entity_pattern = re.compile(r'([A-Za-z0-9\-_,\[\] ]+?)\s*\((human|mouse)\)\s*\(([^\)]+)\)')
        # Explode the first two fields using regex
        for field in [f'Effects of modification on {protein_name}', 'Effects of modification on biological processes']:
            value = downstream_data.get(field, '')
            for match in effect_pattern.finditer(value):
                effect = match.group(1).strip()
                references = [int(ref.strip()) for ref in match.group(2).split(',') if ref.strip().isdigit()]
                downstream_rows.append({
                    'Downstream regulation': field,
                    'Downstream protein': None,
                    'Organism': None,
                    'References': references,
                    'Amino Acid': amino_acid,
                    'Protein': protein_name,
                    'Activity': effect,
                })
        # Explode the last two fields as before
        for field in ['Induce interaction with:', 'Inhibit interaction with:']:
            value = downstream_data.get(field, '')
            for match in entity_pattern.finditer(value):
                protein = match.group(1).strip()
                organism = match.group(2).strip()
                references = [int(ref.strip()) for ref in match.group(3).split(',') if ref.strip().isdigit()]
                downstream_rows.append({
                    'Downstream regulation': field,
                    'Downstream protein': protein.lstrip(', ').strip(),
                    'Organism': organism,
                    'References': references,
                    'Amino Acid': amino_acid,
                    'Protein': protein_name,
                    'Activity': None,
                })

        downstream_columns = [
            'Downstream regulation',
            'Downstream protein',
            'Organism',
            'References',
            'Amino Acid',
            'Protein',
            'Activity',
        ]
        upstream_columns = ['Upstream regulation', 'Upstream protein']
        if not downstream_rows:
            downstream_rows = [{column: None for column in downstream_columns}]

        merged_rows = downstream_rows + upstream_rows
        references = await references_scraper(page)
        pubmed_by_reference = {
            str(reference.get('Reference Number')): reference.get('PubMed ID')
            for reference in references
        }
        output_rows = []
        for row in merged_rows:
            row_references = row.get('References')
            if isinstance(row_references, list) and row_references:
                references_to_write = row_references
            else:
                references_to_write = [None]

            for reference_number in references_to_write:
                output_row = {
                    column: row.get(column)
                    for column in [*downstream_columns, *upstream_columns]
                }
                reference_text = str(reference_number) if reference_number is not None else None
                output_row['Site ID'] = site_id
                output_row['Site URL'] = url
                output_row['Reference Number'] = reference_text
                output_row['PubMed ID'] = pubmed_by_reference.get(reference_text)
                output_rows.append(output_row)

        folder = safe_filename_component(protein_name)
        os.makedirs(folder, exist_ok=True)
        output_filename = os.path.join(
            folder,
            f"{safe_filename_component(amino_acid)}_{safe_filename_component(protein_name)}_site{site_id}.csv",
        )
        fieldnames = [
            'Site ID',
            'Site URL',
            *downstream_columns,
            *upstream_columns,
            'Reference Number',
            'PubMed ID',
        ]
        with open(output_filename, 'w', newline='', encoding='utf-8') as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            for row in output_rows:
                writer.writerow({
                    key: 'nan' if value is None or value == '' else value
                    for key, value in row.items()
                })
        print(f"Saved {output_filename}")
        await close_browser_context(browser, context)


NON_HUMAN_ORGANISM_PATTERN = re.compile(
    r"\b(mouse|rat|pig|cow|bovine|dog|canine|chicken|frog|zebrafish|monkey)\b",
    re.IGNORECASE,
)


def is_human_organism_text(text):
    normalized = " ".join(text.split()).strip().lower()
    return bool(re.search(r"\bhuman\b", normalized)) or "(human)" in normalized


def is_human_site_row(row_text):
    normalized = " ".join(row_text.split())
    return not NON_HUMAN_ORGANISM_PATTERN.search(normalized)


async def collect_human_site_hrefs(page):
    try:
        await page.wait_for_selector("a[href*='siteAction.action?id=']", timeout=30000)
    except Exception:
        pass

    rows = await page.eval_on_selector_all(
        "tr",
        """
        rows => rows.flatMap(row => {
            const rowText = row.textContent || "";
            return Array.from(row.querySelectorAll("a[href*='siteAction.action?id=']")).map(link => ({
                href: link.href,
                rowText
            }));
        })
        """
    )
    hrefs = [row["href"] for row in rows if is_human_site_row(row["rowText"])]
    return hrefs


async def find_site_ids_for_protein(protein_id):
    from playwright.async_api import async_playwright

    protein_url = f"https://www.phosphosite.org/proteinAction.action?id={protein_id}&showAllSites=true"
    site_table_url = f"https://www.phosphosite.org/siteTableNewAction?id={protein_id}&showAllSites=true"
    print(f"STAGE: protein {protein_id}: opening browser for site discovery", flush=True)
    async with async_playwright() as p:
        browser, context, page = await new_browser_context(p)
        try:
            print(f"STAGE: protein {protein_id}: warming PhosphoSitePlus session", flush=True)
            await goto_with_cloudflare_retry(
                page,
                "https://www.phosphosite.org/homeAction",
                "PhosphoSitePlus home",
            )
            print(f"STAGE: protein {protein_id}: loading protein page {protein_url}", flush=True)
            await goto_with_cloudflare_retry(page, protein_url, f"protein ID {protein_id}")
            hrefs = await collect_human_site_hrefs(page)
            print(f"STAGE: protein {protein_id}: found {len(hrefs)} human site link(s) on protein page", flush=True)

            if not hrefs:
                print(f"STAGE: protein {protein_id}: trying site table fallback {site_table_url}", flush=True)
                await goto_with_cloudflare_retry(page, site_table_url, f"site table protein ID {protein_id}")
                hrefs = await collect_human_site_hrefs(page)
                print(f"STAGE: protein {protein_id}: found {len(hrefs)} human site link(s) in site table", flush=True)
        finally:
            print(f"STAGE: protein {protein_id}: closing site discovery browser", flush=True)
            await close_browser_context(browser, context)

    site_ids = []
    for href in hrefs:
        match = re.search(r"siteAction\.action\?id=(\d+)", href)
        if match:
            site_ids.append(int(match.group(1)))
    return list(dict.fromkeys(site_ids))


async def scrape_protein(protein_id, delay, continue_on_error, delay_jitter=DEFAULT_DELAY_JITTER):
    print(f"START: protein ID {protein_id}: discovering human siteAction IDs", flush=True)
    site_ids = await find_site_ids_for_protein(protein_id)
    if not site_ids:
        raise RuntimeError(f"No siteAction IDs found for protein ID {protein_id}.")

    print(f"FOUND: protein ID {protein_id}: {len(site_ids)} unique human site(s)", flush=True)
    print(f"START: protein ID {protein_id}: scraping site IDs", flush=True)
    site_result = await run_site_batch(site_ids, delay, continue_on_error, delay_jitter)
    print(f"DONE: protein ID {protein_id}: completed site scrape batch", flush=True)
    return site_result


def extract_protein_id_from_url(url):
    match = re.search(r"proteinAction\.action\?id=(\d+)", url)
    if match:
        return int(match.group(1))
    return None


async def get_selected_protein_organism(page):
    try:
        selected = await page.eval_on_selector(
            "#select_id",
            """
            select => {
                if (!select || select.selectedIndex < 0) return "";
                return select.options[select.selectedIndex].textContent.trim();
            }
            """,
        )
        return selected or ""
    except Exception:
        return ""


async def navigate_to_protein_url(page, url, require_human=True):
    print(f"STAGE: navigating to candidate protein URL {url}", flush=True)
    await goto_with_cloudflare_retry(page, url, "candidate protein URL")
    protein_id = extract_protein_id_from_url(page.url)
    if protein_id is None:
        raise RuntimeError(f"Navigation did not land on a proteinAction URL: {page.url}")
    print(f"STAGE: landed on final protein URL {page.url}", flush=True)
    if require_human:
        organism = await get_selected_protein_organism(page)
        if organism and not is_human_organism_text(organism):
            raise RuntimeError(f"Resolved URL is not human: organism='{organism}', url={page.url}")
        print(f"STAGE: final URL organism check passed ({organism or 'human assumed from result row'})", flush=True)
    return protein_id, page.url


def result_from_human_search_candidate(protein_name, protein_id, url, source):
    print(
        f"RESOLVED: {protein_name}: id={protein_id} from human search-result URL {url}",
        flush=True,
    )
    return {
        "protein_name": protein_name,
        "protein_id": protein_id,
        "url": url,
        "organism": "human",
        "source": source,
    }


async def resolve_protein_name_on_page(page, protein_name, organism):
    if organism.strip().lower() != "human":
        raise ValueError("This workflow only supports human protein lookup.")

    search_url = (
        "https://www.phosphosite.org/simpleSearchSubmitAction.action"
        f"?searchStr={quote_plus(protein_name)}"
    )
    print(f"STAGE: {protein_name}: warming PhosphoSitePlus home page", flush=True)
    await goto_with_cloudflare_retry(
        page,
        "https://www.phosphosite.org/homeAction",
        f"{protein_name} home warmup",
    )
    print(f"STAGE: {protein_name}: searching {search_url}", flush=True)
    await goto_with_cloudflare_retry(page, search_url, f"{protein_name} search")

    protein_id = extract_protein_id_from_url(page.url)
    if protein_id is not None:
        print(f"STAGE: {protein_name}: search landed directly on protein page {page.url}", flush=True)
        selected_organism = await get_selected_protein_organism(page)
        if selected_organism and not is_human_organism_text(selected_organism):
            raise RuntimeError(
                f"Direct lookup for {protein_name} resolved to non-human organism '{selected_organism}'."
            )
        print(f"RESOLVED: {protein_name}: id={protein_id} from final URL {page.url}", flush=True)
        return {
            "protein_name": protein_name,
            "protein_id": protein_id,
            "url": page.url,
            "organism": "human",
            "source": "direct_url",
        }

    try:
        await page.wait_for_selector(
            "#simpleSearchResultsTable tbody tr a[href*='proteinAction.action?id=']",
            timeout=20000,
        )
    except Exception:
        pass

    rows = await page.eval_on_selector_all(
        "#simpleSearchResultsTable tbody tr",
        """
        rows => rows.map(row => {
            const cells = Array.from(row.querySelectorAll("td")).map(cell => cell.textContent.trim());
            const links = Array.from(row.querySelectorAll("a[href*='proteinAction.action?id=']")).map(link => ({
                href: link.href,
                text: link.textContent.trim()
            }));
            return {
                protein: cells[0] || "",
                gene: cells[1] || "",
                organism: cells[2] || "",
                links
            };
        })
        """
    )

    normalized_name = protein_name.strip().lower()
    candidates = []
    human_result_rows = 0
    for row in rows:
        if is_human_organism_text(row["organism"].strip().lower()):
            human_result_rows += 1
        for link in row["links"]:
            match = re.search(r"proteinAction\.action\?id=(\d+)", link["href"])
            if not match:
                continue

            gene_tokens = [
                token.strip().lower()
                for token in re.split(r"[,;/\s]+", row["gene"])
                if token.strip()
            ]
            protein_text = row["protein"].strip().lower()
            link_text = link["text"].strip().lower()
            organism_text = row["organism"].strip().lower()
            if not is_human_organism_text(organism_text):
                continue

            score = 0
            if normalized_name in gene_tokens:
                score += 100
            if protein_text == normalized_name:
                score += 80
            if link_text == normalized_name:
                score += 80
            if normalized_name in protein_text or normalized_name in link_text:
                score += 10

            candidates.append((score, link["href"], row))

    print(
        f"STAGE: {protein_name}: inspected {len(rows)} search result row(s), "
        f"{human_result_rows} human row(s), {len(candidates)} human protein URL candidate(s)",
        flush=True,
    )
    if candidates:
        candidates.sort(key=lambda candidate: candidate[0], reverse=True)
        best_score, best_url, best_row = candidates[0]
        if best_score > 0:
            print(
                f"STAGE: {protein_name}: selecting best human candidate score={best_score}, "
                f"gene='{best_row['gene']}', protein='{best_row['protein']}'",
                flush=True,
            )
            best_id_from_url = extract_protein_id_from_url(best_url)
            try:
                best_id, resolved_url = await navigate_to_protein_url(page, best_url)
            except RuntimeError as exc:
                final_url = getattr(page, "url", "")
                if best_id_from_url is not None and "confirmAction" in final_url:
                    print(
                        f"WARNING: {protein_name}: candidate URL redirected to confirmAction; "
                        f"saving human search-result id={best_id_from_url} from {best_url}",
                        flush=True,
                    )
                    return result_from_human_search_candidate(
                        protein_name,
                        best_id_from_url,
                        best_url,
                        "human_search_result_url_confirm_redirect",
                    )
                raise exc
            print(
                f"RESOLVED: {protein_name}: id={best_id} "
                f"from URL: {resolved_url}; search row: protein='{best_row['protein']}', "
                f"gene='{best_row['gene']}', organism='{best_row['organism']}'",
                flush=True,
            )
            return {
                "protein_name": protein_name,
                "protein_id": best_id,
                "url": resolved_url,
                "organism": "human",
                "source": "human_search_result_final_url",
            }

    raise RuntimeError(f"No human proteinAction ID found for protein name '{protein_name}'.")


async def resolve_protein_name(protein_name, organism):
    from playwright.async_api import async_playwright

    async with async_playwright() as p:
        print(f"STAGE: {protein_name}: launching browser", flush=True)
        browser, context, page = await new_browser_context(p)
        try:
            return await resolve_protein_name_on_page(page, protein_name, organism)
        finally:
            print(f"STAGE: {protein_name}: closing browser", flush=True)
            await close_browser_context(browser, context)


async def resolve_protein_names(protein_names, organism, delay, continue_on_error=False, delay_jitter=DEFAULT_DELAY_JITTER):
    from playwright.async_api import async_playwright

    resolved = []
    failures = []

    async with async_playwright() as p:
        print("SESSION: launching one browser context for lookup batch", flush=True)
        browser, context, page = await new_browser_context(p)
        try:
            for index, protein_name in enumerate(protein_names, start=1):
                print(f"[{index}/{len(protein_names)}] START lookup {protein_name}", flush=True)
                try:
                    result = await resolve_protein_name_on_page(page, protein_name, organism)
                    resolved.append(result)
                    print(
                        f"[{index}/{len(protein_names)}] DONE lookup {protein_name}: "
                        f"id={result['protein_id']} url={result['url']}",
                        flush=True,
                    )
                except Exception as exc:
                    failures.append((protein_name, exc))
                    print(f"ERROR: {protein_name} lookup failed: {exc}", flush=True)
                    if not continue_on_error:
                        raise

                if index < len(protein_names):
                    await sleep_between_requests(delay, delay_jitter, "protein lookup")
        finally:
            print("SESSION: closing lookup browser context", flush=True)
            await close_browser_context(browser, context)

    if failures:
        failed_names = ", ".join(name for name, _ in failures)
        print(f"FAILED: protein name lookups: {failed_names}", flush=True)

    return resolved


def write_lookup_results(path, results):
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["protein_name", "protein_id", "url", "organism", "source"],
        )
        writer.writeheader()
        writer.writerows(results)


def parse_ids_file(path, label):
    """Read numeric IDs from a curator-maintained text file."""
    ids = []
    with open(path, "r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            value = line.split("#", 1)[0].strip()
            if not value:
                continue
            try:
                ids.append(int(value))
            except ValueError as exc:
                raise ValueError(f"Invalid {label} on line {line_number}: {value}") from exc
    return ids


def parse_names_file(path):
    """Read protein names from a curator-maintained text file."""
    names = []
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            value = line.split("#", 1)[0].strip().rstrip(",").strip("'\" \t")
            if value:
                names.append(value)
    return names


def build_parser():
    parser = argparse.ArgumentParser(
        description="Scrape PhosphoSitePlus phosphorylation site data."
    )
    parser.add_argument(
        "protein_id",
        nargs="?",
        type=int,
        help="Single PhosphoSitePlus protein ID to scrape.",
    )
    parser.add_argument(
        "--protein-id",
        dest="protein_id_option",
        type=int,
        help="Single PhosphoSitePlus protein ID to scrape.",
    )
    parser.add_argument(
        "--protein-ids",
        nargs="+",
        type=int,
        help="One or more PhosphoSitePlus protein IDs to scrape.",
    )
    parser.add_argument(
        "--protein-ids-file",
        help="Text file containing one protein ID per line. Comments after # are ignored.",
    )
    parser.add_argument(
        "--protein-name",
        dest="protein_name_option",
        help="Single protein name/symbol to look up and scrape.",
    )
    parser.add_argument(
        "--protein-names",
        nargs="+",
        help="One or more protein names/symbols to look up and scrape.",
    )
    parser.add_argument(
        "--protein-names-file",
        help="Text file containing one protein name/symbol per line.",
    )
    parser.add_argument(
        "--organism",
        default="human",
        help="Protein lookup organism. This workflow only supports human.",
    )
    parser.add_argument(
        "--lookup-only",
        action="store_true",
        help="Resolve protein names to IDs and exit without scraping.",
    )
    parser.add_argument(
        "--lookup-output",
        help="Optional CSV path for resolved protein_name, protein_id, url, and source.",
    )
    parser.add_argument(
        "--site-id",
        dest="site_id_option",
        type=int,
        help="Advanced: scrape a single siteAction SITE_ID directly.",
    )
    parser.add_argument(
        "--site-ids",
        nargs="+",
        type=int,
        help="Advanced: scrape one or more siteAction SITE_IDs directly.",
    )
    parser.add_argument(
        "--site-ids-file",
        help="Advanced: text file containing one siteAction SITE_ID per line.",
    )
    parser.add_argument(
        "--delay",
        type=float,
        default=2.0,
        help="Seconds to wait between IDs in batch mode. Default: 2.0.",
    )
    parser.add_argument(
        "--delay-jitter",
        type=float,
        default=DEFAULT_DELAY_JITTER,
        help="Randomize waits by this fraction around --delay. Example: 0.5 gives 50%%-150%% of delay.",
    )
    parser.add_argument(
        "--continue-on-error",
        action="store_true",
        help="Continue a batch if one ID fails.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Validate and print input IDs without scraping.",
    )
    return parser


async def run_site_batch(site_ids, delay, continue_on_error, delay_jitter=DEFAULT_DELAY_JITTER):
    if not site_ids:
        raise ValueError("No SITE_IDs were provided.")

    print(f"START: scraping {len(site_ids)} siteAction ID(s)", flush=True)
    failures = []
    for index, site_id in enumerate(site_ids, start=1):
        print(f"[{index}/{len(site_ids)}] START siteAction ID {site_id}", flush=True)
        try:
            await main(site_id)
            print(f"[{index}/{len(site_ids)}] DONE siteAction ID {site_id}", flush=True)
        except Exception as exc:
            failures.append((site_id, exc))
            print(f"ERROR: SITE_ID {site_id} failed: {exc}", flush=True)
            if not continue_on_error:
                break

        if index < len(site_ids):
            await sleep_between_requests(delay, delay_jitter, "siteAction scrape")

    if failures:
        failed_ids = ", ".join(str(site_id) for site_id, _ in failures)
        if continue_on_error:
            print(f"WARNING: failed SITE_IDs after retries: {failed_ids}", flush=True)
            return {
                "failed_site_ids": [site_id for site_id, _ in failures],
                "status": "completed_with_site_errors",
            }
        raise RuntimeError(f"Failed SITE_IDs: {failed_ids}")
    print("DONE: siteAction scrape batch completed", flush=True)
    return {"failed_site_ids": [], "status": "complete"}


async def run_protein_batch(protein_ids, delay, continue_on_error, delay_jitter=DEFAULT_DELAY_JITTER):
    if not protein_ids:
        raise ValueError("No protein IDs were provided.")

    failures = []
    site_failures_by_protein = {}
    for index, protein_id in enumerate(protein_ids, start=1):
        print(f"[{index}/{len(protein_ids)}] START protein ID {protein_id}", flush=True)
        try:
            site_result = await scrape_protein(protein_id, delay, continue_on_error, delay_jitter)
            failed_site_ids = site_result.get("failed_site_ids", []) if site_result else []
            if failed_site_ids:
                site_failures_by_protein[str(protein_id)] = failed_site_ids
            print(f"[{index}/{len(protein_ids)}] DONE protein ID {protein_id}", flush=True)
        except Exception as exc:
            failures.append((protein_id, exc))
            print(f"ERROR: Protein ID {protein_id} failed: {exc}", flush=True)
            if not continue_on_error:
                break

        if index < len(protein_ids):
            await sleep_between_requests(delay, delay_jitter, "protein ID")

    if failures:
        failed_ids = ", ".join(str(protein_id) for protein_id, _ in failures)
        if continue_on_error:
            print(f"WARNING: failed protein IDs after retries: {failed_ids}", flush=True)
            return {
                "failed_protein_ids": [protein_id for protein_id, _ in failures],
                "failed_site_ids_by_protein": site_failures_by_protein,
                "status": "completed_with_protein_errors",
            }
        raise RuntimeError(f"Failed protein IDs: {failed_ids}")
    print("DONE: protein scrape batch completed", flush=True)
    status = "completed_with_site_errors" if site_failures_by_protein else "complete"
    return {
        "failed_protein_ids": [],
        "failed_site_ids_by_protein": site_failures_by_protein,
        "status": status,
    }


def collect_protein_ids(args):
    protein_ids = []
    if args.protein_id is not None:
        protein_ids.append(args.protein_id)
    if args.protein_id_option is not None:
        protein_ids.append(args.protein_id_option)
    if args.protein_ids:
        protein_ids.extend(args.protein_ids)
    if args.protein_ids_file:
        protein_ids.extend(parse_ids_file(args.protein_ids_file, "protein ID"))
    return protein_ids


def collect_protein_names(args):
    protein_names = []
    if args.protein_name_option:
        protein_names.append(args.protein_name_option)
    if args.protein_names:
        protein_names.extend(args.protein_names)
    if args.protein_names_file:
        protein_names.extend(parse_names_file(args.protein_names_file))
    return protein_names


def collect_site_ids(args):
    site_ids = []
    if args.site_id_option is not None:
        site_ids.append(args.site_id_option)
    if args.site_ids:
        site_ids.extend(args.site_ids)
    if args.site_ids_file:
        site_ids.extend(parse_ids_file(args.site_ids_file, "SITE_ID"))
    return site_ids


if __name__ == "__main__":
    args = build_parser().parse_args()
    protein_ids = collect_protein_ids(args)
    protein_names = collect_protein_names(args)
    site_ids = collect_site_ids(args)

    if not protein_ids and not protein_names and not site_ids:
        protein_input = input("Enter protein name or proteinAction ID (e.g., CDK1 or 582): ")
        try:
            protein_ids = [int(protein_input)]
        except ValueError:
            protein_names = [protein_input]

    if args.dry_run:
        if protein_names:
            print(f"Found {len(protein_names)} protein name(s): {', '.join(protein_names)}")
        if protein_ids:
            print(f"Found {len(protein_ids)} protein ID(s): {', '.join(str(protein_id) for protein_id in protein_ids)}")
        if site_ids:
            print(f"Found {len(site_ids)} direct SITE_ID(s): {', '.join(str(site_id) for site_id in site_ids)}")
        exit(0)

    if protein_names:
        resolved = asyncio.run(
            resolve_protein_names(
                protein_names,
                args.organism,
                args.delay,
                args.continue_on_error,
                args.delay_jitter,
            )
        )
        if args.lookup_output:
            write_lookup_results(args.lookup_output, resolved)
            print(f"Saved lookup results to {args.lookup_output}")
        protein_ids.extend(result["protein_id"] for result in resolved)
        if args.lookup_only:
            exit(0)

    if protein_ids:
        asyncio.run(run_protein_batch(protein_ids, args.delay, args.continue_on_error, args.delay_jitter))
    if site_ids:
        asyncio.run(run_site_batch(site_ids, args.delay, args.continue_on_error, args.delay_jitter))
