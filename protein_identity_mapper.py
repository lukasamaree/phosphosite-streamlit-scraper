import argparse
import csv
import json
import re
import ssl
import time
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import quote, urlencode
from urllib.request import Request, urlopen


ROOT = Path(__file__).resolve().parent
CURATION_DIR = ROOT / "curated_protein_ids"
DEFAULT_OUTPUT = CURATION_DIR / "protein_identity_lookup.csv"
DEFAULT_CACHE = CURATION_DIR / "protein_identity_lookup_cache.json"
DEFAULT_COLUMNS = ("Protein", "Downstream protein", "Upstream protein")
HUMAN_ORGANISM_ID = "9606"


def normalize_name(value):
    text = str(value or "").strip()
    text = re.sub(r"\([^)]*\b(human|mouse|rat|rabbit|pig|hamster|fruit fly)[^)]*\)", "", text, flags=re.I)
    text = re.sub(r"\biso\d+\b", "", text, flags=re.I)
    text = re.sub(r"\s+", " ", text.replace("\xa0", " ")).strip(" ,;")
    return text


def name_key(value):
    return re.sub(r"[^a-z0-9]+", "", normalize_name(value).lower())


def unique_ordered(values):
    seen = set()
    output = []
    for value in values:
        normalized = normalize_name(value)
        key = name_key(normalized)
        if normalized and key and key not in seen:
            seen.add(key)
            output.append(normalized)
    return output


def split_aliases(value):
    if isinstance(value, list):
        items = value
    else:
        items = re.split(r"[|;,]", str(value or ""))
    return unique_ordered(items)


def request_json(url, headers=None, timeout=30):
    request = Request(url, headers=headers or {"User-Agent": "phosphosite-identity-mapper/1.0"})
    try:
        with urlopen(request, timeout=timeout) as response:
            return json.loads(response.read().decode("utf-8"))
    except URLError as exc:
        if "CERTIFICATE_VERIFY_FAILED" not in str(exc):
            raise
        context = ssl._create_unverified_context()
        with urlopen(request, timeout=timeout, context=context) as response:
            return json.loads(response.read().decode("utf-8"))


def safe_request_json(url, headers=None):
    try:
        return request_json(url, headers=headers)
    except (HTTPError, URLError, TimeoutError, json.JSONDecodeError) as exc:
        return {"_error": str(exc)}


def uniprot_query(query, size=5):
    params = {
        "query": f'(organism_id:{HUMAN_ORGANISM_ID}) AND (gene_exact:"{query}" OR gene:"{query}" OR protein_name:"{query}")',
        "format": "json",
        "size": str(size),
    }
    url = "https://rest.uniprot.org/uniprotkb/search?" + urlencode(params)
    return safe_request_json(url)


def hgnc_query(query):
    escaped = quote(query)
    headers = {
        "Accept": "application/json",
        "User-Agent": "phosphosite-identity-mapper/1.0",
    }
    fields = ("symbol", "alias_symbol", "prev_symbol", "name", "alias_name", "prev_name")
    docs = []
    for field in fields:
        url = f"https://rest.genenames.org/search/{field}/{escaped}"
        payload = safe_request_json(url, headers=headers)
        docs.extend(payload.get("response", {}).get("docs", []))
    return docs


def ncbi_gene_query(query):
    search_params = {
        "db": "gene",
        "term": f"{query}[All Fields] AND human[Organism]",
        "retmode": "json",
        "retmax": "5",
    }
    search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?" + urlencode(search_params)
    search = safe_request_json(search_url)
    ids = search.get("esearchresult", {}).get("idlist", [])
    if not ids:
        return []

    summary_params = {"db": "gene", "id": ",".join(ids), "retmode": "json"}
    summary_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?" + urlencode(summary_params)
    summary = safe_request_json(summary_url)
    result = summary.get("result", {})
    return [result[gene_id] for gene_id in ids if gene_id in result]


def ensembl_xref_query(query):
    url = f"https://rest.ensembl.org/xrefs/symbol/homo_sapiens/{quote(query)}?content-type=application/json"
    payload = safe_request_json(url, headers={"User-Agent": "phosphosite-identity-mapper/1.0"})
    return payload if isinstance(payload, list) else []


def extract_uniprot_hit(hit):
    genes = hit.get("genes", [])
    primary_gene = ""
    gene_synonyms = []
    if genes:
        primary_gene = genes[0].get("geneName", {}).get("value", "")
        for gene in genes:
            gene_synonyms.extend(item.get("value", "") for item in gene.get("synonyms", []))

    protein = hit.get("proteinDescription", {})
    recommended_name = protein.get("recommendedName", {}).get("fullName", {}).get("value", "")
    alternative_names = [
        item.get("fullName", {}).get("value", "")
        for item in protein.get("alternativeNames", [])
    ]
    return {
        "canonical_gene": primary_gene,
        "uniprot_accession": hit.get("primaryAccession", ""),
        "recommended_name": recommended_name,
        "aliases": unique_ordered([primary_gene, *gene_synonyms, recommended_name, *alternative_names]),
    }


def best_hgnc_doc(query, docs):
    query_key = name_key(query)
    for doc in docs:
        candidates = [doc.get("symbol", ""), *doc.get("alias_symbol", []), *doc.get("prev_symbol", [])]
        if query_key in {name_key(candidate) for candidate in candidates}:
            return doc
    return docs[0] if docs else {}


def best_ncbi_doc(query, docs):
    query_key = name_key(query)
    for doc in docs:
        candidates = [
            doc.get("name", ""),
            doc.get("nomenclaturesymbol", ""),
            *split_aliases(doc.get("otheraliases", "")),
        ]
        if query_key in {name_key(candidate) for candidate in candidates}:
            return doc
    return {}


def source_match(query, aliases):
    query_key = name_key(query)
    return query_key in {name_key(alias) for alias in aliases}


def confidence_for(query, canonical_gene, uniprot_accession, aliases, sources):
    if not canonical_gene:
        return "unresolved"
    if uniprot_accession and source_match(query, aliases):
        return "high"
    if uniprot_accession:
        return "medium"
    if "HGNC" in sources or "NCBI" in sources:
        return "medium"
    return "low"


def resolve_identity(raw_name, cache=None, sleep_seconds=0.1):
    normalized_query = normalize_name(raw_name)
    cache = cache if cache is not None else {}
    cache_key = normalized_query.lower()
    if cache_key in cache:
        return cache[cache_key]

    uniprot_payload = uniprot_query(normalized_query)
    hgnc_docs = hgnc_query(normalized_query)
    ncbi_docs = ncbi_gene_query(normalized_query)
    ensembl_docs = ensembl_xref_query(normalized_query)

    sources = []
    canonical_gene = ""
    uniprot_accession = ""
    recommended_name = ""
    aliases = [normalized_query]
    hgnc_id = ""
    ncbi_gene_id = ""
    ensembl_gene_id = ""

    hgnc_doc = best_hgnc_doc(normalized_query, hgnc_docs)
    if hgnc_doc:
        sources.append("HGNC")
        canonical_gene = hgnc_doc.get("symbol", "") or canonical_gene
        hgnc_id = hgnc_doc.get("hgnc_id", "")
        recommended_name = hgnc_doc.get("name", "") or recommended_name
        uniprot_ids = hgnc_doc.get("uniprot_ids", [])
        if uniprot_ids:
            uniprot_accession = uniprot_ids[0]
        aliases.extend(
            [
                canonical_gene,
                hgnc_doc.get("name", ""),
                *hgnc_doc.get("alias_symbol", []),
                *hgnc_doc.get("prev_symbol", []),
                *hgnc_doc.get("alias_name", []),
                *hgnc_doc.get("prev_name", []),
            ]
        )

    uniprot_results = uniprot_payload.get("results", [])
    if uniprot_results:
        sources.append("UniProt")
        uniprot_hit = extract_uniprot_hit(uniprot_results[0])
        canonical_gene = canonical_gene or uniprot_hit["canonical_gene"]
        uniprot_accession = uniprot_accession or uniprot_hit["uniprot_accession"]
        recommended_name = recommended_name or uniprot_hit["recommended_name"]
        aliases.extend(uniprot_hit["aliases"])

    if canonical_gene and not uniprot_accession:
        canonical_payload = uniprot_query(canonical_gene, size=1)
        canonical_results = canonical_payload.get("results", [])
        if canonical_results:
            sources.append("UniProt")
            canonical_hit = extract_uniprot_hit(canonical_results[0])
            uniprot_accession = canonical_hit["uniprot_accession"]
            recommended_name = recommended_name or canonical_hit["recommended_name"]
            aliases.extend(canonical_hit["aliases"])

    ncbi_doc = best_ncbi_doc(normalized_query, ncbi_docs)
    if ncbi_doc:
        sources.append("NCBI")
        canonical_gene = canonical_gene or ncbi_doc.get("nomenclaturesymbol") or ncbi_doc.get("name", "")
        recommended_name = recommended_name or ncbi_doc.get("description", "")
        ncbi_gene_id = str(ncbi_doc.get("uid", ""))
        aliases.extend([ncbi_doc.get("name", ""), ncbi_doc.get("description", ""), *split_aliases(ncbi_doc.get("otheraliases", ""))])

    if ensembl_docs:
        sources.append("Ensembl")
        ensembl_gene = next((doc for doc in ensembl_docs if doc.get("type") == "gene"), ensembl_docs[0])
        ensembl_gene_id = ensembl_gene.get("id", "")
        aliases.extend([ensembl_gene.get("display_id", ""), ensembl_gene.get("description", "")])

    aliases = unique_ordered(aliases)
    sources = unique_ordered(sources)
    confidence = confidence_for(normalized_query, canonical_gene, uniprot_accession, aliases, sources)
    result = {
        "raw_name": raw_name,
        "normalized_query": normalized_query,
        "canonical_gene": canonical_gene,
        "uniprot_accession": uniprot_accession,
        "recommended_name": recommended_name,
        "organism": "human",
        "aliases": ";".join(aliases),
        "confidence": confidence,
        "match_status": "resolved" if canonical_gene or uniprot_accession else "unresolved",
        "sources": ";".join(sources),
        "hgnc_id": hgnc_id,
        "ncbi_gene_id": ncbi_gene_id,
        "ensembl_gene_id": ensembl_gene_id,
    }
    cache[cache_key] = result
    if sleep_seconds:
        time.sleep(sleep_seconds)
    return result


def read_names_file(path):
    names = []
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            value = line.split("#", 1)[0].strip().strip(",")
            if value:
                names.append(value.strip("'\""))
    return names


def collect_names_from_outputs(root, columns=DEFAULT_COLUMNS):
    names = []
    columns_lower = {column.lower(): column for column in columns}
    for path in Path(root).rglob("*.csv"):
        if any(part.startswith(".") or part == "__pycache__" for part in path.parts):
            continue
        try:
            with open(path, "r", encoding="utf-8", errors="replace", newline="") as handle:
                reader = csv.DictReader(handle)
                field_map = {field.lower(): field for field in reader.fieldnames or []}
                selected_fields = [field_map[key] for key in columns_lower if key in field_map]
                for row in reader:
                    for field in selected_fields:
                        value = normalize_name(row.get(field))
                        if value:
                            names.append(value)
        except Exception:
            continue
    return unique_ordered(names)


def load_cache(path, refresh=False):
    path = Path(path)
    if refresh or not path.exists():
        return {}
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def save_cache(path, cache):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(cache, handle, indent=2, sort_keys=True)


def write_rows(path, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
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
    ]
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def build_parser():
    parser = argparse.ArgumentParser(description="Build a canonical gene/UniProt lookup table for PhosphoSitePlus protein names.")
    parser.add_argument("--names", nargs="+", help="Raw protein names to resolve.")
    parser.add_argument("--names-file", help="Text file with one raw protein name per line.")
    parser.add_argument("--output-root", default=str(ROOT), help="Scraper output root to scan for CSV protein columns.")
    parser.add_argument("--columns", nargs="+", default=list(DEFAULT_COLUMNS), help="CSV columns to scan.")
    parser.add_argument("--output-csv", default=str(DEFAULT_OUTPUT), help="Lookup CSV to write.")
    parser.add_argument("--cache-json", default=str(DEFAULT_CACHE), help="API result cache JSON.")
    parser.add_argument("--refresh-cache", action="store_true", help="Ignore cached resolved names.")
    return parser


def main():
    args = build_parser().parse_args()
    names = []
    if args.names:
        names.extend(args.names)
    if args.names_file:
        names.extend(read_names_file(args.names_file))
    if not names:
        names.extend(collect_names_from_outputs(args.output_root, args.columns))
    names = unique_ordered(names)
    if not names:
        raise SystemExit("No protein names found.")

    cache = load_cache(args.cache_json, args.refresh_cache)
    rows = [resolve_identity(name, cache=cache) for name in names]
    write_rows(args.output_csv, rows)
    save_cache(args.cache_json, cache)

    resolved = sum(1 for row in rows if row["match_status"] == "resolved")
    print(f"IDENTITY_LOOKUP: wrote {len(rows)} row(s) to {args.output_csv}", flush=True)
    print(f"IDENTITY_LOOKUP: resolved={resolved}, unresolved={len(rows) - resolved}", flush=True)
    print("IDENTITY_LOOKUP_JSON: " + json.dumps({"rows": len(rows), "resolved": resolved, "output_csv": str(Path(args.output_csv).resolve())}, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
