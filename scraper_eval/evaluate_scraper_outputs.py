import argparse
import csv
import json
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path


REQUIRED_OUTPUT_COLUMNS = {
    "Protein",
    "Amino Acid",
    "Organism",
    "References",
    "PubMed ID",
}


@dataclass(frozen=True)
class ExpectedProtein:
    query: str
    expected_protein_id: str
    expected_organism: str
    expected_uniprot: str
    output_dir: str
    required_sites: tuple[str, ...]


def normalize(value):
    return str(value or "").strip()


def normalize_key(value):
    return normalize(value).lower()


def is_blankish(value):
    return normalize_key(value) in {"", "nan", "none", "null"}


def read_csv_rows(path):
    with open(path, "r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle)
        return reader.fieldnames or [], list(reader)


def load_expected_manifest(path):
    _, rows = read_csv_rows(path)
    expected = []
    for row in rows:
        required_sites = tuple(
            site.strip()
            for site in normalize(row.get("required_sites")).split(";")
            if site.strip()
        )
        expected.append(
            ExpectedProtein(
                query=normalize(row.get("query")),
                expected_protein_id=normalize(row.get("expected_protein_id")),
                expected_organism=normalize(row.get("expected_organism")).lower(),
                expected_uniprot=normalize(row.get("expected_uniprot")),
                output_dir=normalize(row.get("output_dir") or row.get("query")),
                required_sites=required_sites,
            )
        )
    return expected


def load_curated_ids(path):
    if not path.exists():
        return {}
    _, rows = read_csv_rows(path)
    by_name = {}
    for row in rows:
        protein_name = normalize(row.get("protein_name"))
        if protein_name:
            by_name[normalize_key(protein_name)] = row
    return by_name


def find_output_files(output_root, output_dir):
    target_dir = output_root / output_dir
    if not target_dir.exists():
        return []
    return sorted(target_dir.glob("*.csv"), key=lambda item: item.name.lower())


def collect_output_rows(files):
    output_rows = []
    schema_errors = []
    for path in files:
        try:
            columns, rows = read_csv_rows(path)
        except Exception as exc:
            schema_errors.append(f"{path}: could not read CSV: {exc}")
            continue

        missing = sorted(REQUIRED_OUTPUT_COLUMNS - set(columns))
        if missing:
            schema_errors.append(f"{path}: missing required columns: {', '.join(missing)}")

        for row in rows:
            row["_source_file"] = str(path)
            output_rows.append(row)
    return output_rows, schema_errors


def site_key(row):
    return (
        normalize_key(row.get("Protein")),
        normalize_key(row.get("Amino Acid")),
        normalize_key(row.get("Organism")),
        normalize(row.get("PubMed ID")),
    )


def extract_sites(rows):
    return {
        normalize(row.get("Amino Acid"))
        for row in rows
        if not is_blankish(row.get("Amino Acid"))
    }


def is_placeholder_ptm_row(row):
    return all(is_blankish(row.get(column)) for column in REQUIRED_OUTPUT_COLUMNS)


def duplicate_rows(rows):
    counts = Counter(site_key(row) for row in rows if not is_placeholder_ptm_row(row))
    return [key for key, count in counts.items() if count > 1]


def classify_result(result):
    signals = []
    action = "use_output"

    if result["critical_identity_failure"]:
        signals.append("critical_identity_failure")
        return "IDENTITY_FAILURE", signals, "stop_and_resolve_identity"

    if result["schema_errors"]:
        signals.append("schema_errors")
    if result["output_row_count"] == 0:
        signals.append("empty_output")
    if result["placeholder_row_count"] > 0:
        signals.append("placeholder_rows_detected")
    if result["missing_required_sites"]:
        signals.append("missing_required_sites")
    if result["duplicate_row_count"] > 0:
        signals.append("duplicate_rows")

    if "schema_errors" in signals:
        failure_class = "SCHEMA_FAILURE"
        action = "inspect_parser_or_output_shape"
    elif "empty_output" in signals:
        failure_class = "EMPTY_OUTPUT"
        action = "retry_later_with_delay"
    elif "placeholder_rows_detected" in signals:
        failure_class = "CLOUDFLARE_LIKELY"
        action = "retry_later_with_delay"
    elif "missing_required_sites" in signals or "duplicate_rows" in signals:
        failure_class = "PARTIAL_SCRAPE"
        action = "retry_missing_or_review_output"
    else:
        failure_class = "VALID_SCRAPE"

    return failure_class, signals, action


def identity_result(expected, curated_row, output_rows):
    errors = []
    warnings = []

    if not curated_row:
        errors.append(f"{expected.query}: missing curated ID row")
        return {
            "identity_score": 0,
            "critical_identity_failure": True,
            "errors": errors,
            "warnings": warnings,
        }

    actual_id = normalize(curated_row.get("protein_id"))
    actual_organism = normalize(curated_row.get("organism")).lower()
    actual_uniprot = normalize(curated_row.get("uniprot") or curated_row.get("accession") or curated_row.get("swissprot"))

    if actual_id != expected.expected_protein_id:
        errors.append(
            f"{expected.query}: wrong protein_id expected={expected.expected_protein_id} actual={actual_id}"
        )
    if actual_organism != expected.expected_organism:
        errors.append(
            f"{expected.query}: wrong organism expected={expected.expected_organism} actual={actual_organism}"
        )
    if expected.expected_uniprot and actual_uniprot != expected.expected_uniprot:
        errors.append(
            f"{expected.query}: wrong UniProt expected={expected.expected_uniprot} actual={actual_uniprot or '<missing>'}"
        )

    output_organisms = {
        normalize_key(row.get("Organism"))
        for row in output_rows
        if not is_blankish(row.get("Organism"))
    }
    non_human = sorted(value for value in output_organisms if value and value != expected.expected_organism)
    if non_human:
        errors.append(f"{expected.query}: output contains non-target organism values: {non_human}")

    return {
        "identity_score": 0 if errors else 1,
        "critical_identity_failure": bool(errors),
        "errors": errors,
        "warnings": warnings,
    }


def ptm_result(expected, output_rows, schema_errors):
    errors = list(schema_errors)
    warnings = []
    observed_sites = extract_sites(output_rows)
    missing_sites = [site for site in expected.required_sites if site not in observed_sites]
    placeholder_rows = [row for row in output_rows if is_placeholder_ptm_row(row)]
    duplicates = duplicate_rows(output_rows)

    if not output_rows:
        errors.append(f"{expected.query}: no PTM output rows found")
    if placeholder_rows:
        errors.append(f"{expected.query}: placeholder PTM rows detected: {len(placeholder_rows)}")
    if missing_sites:
        errors.append(f"{expected.query}: missing required PTM sites: {', '.join(missing_sites)}")
    if duplicates:
        errors.append(f"{expected.query}: duplicate PTM rows detected: {len(duplicates)}")

    checks = 4
    failures = (
        int(bool(schema_errors))
        + int(bool(placeholder_rows or not output_rows))
        + int(bool(missing_sites or not output_rows))
        + int(bool(duplicates))
    )
    ptm_score = max(0.0, (checks - failures) / checks)
    return {
        "ptm_score": ptm_score,
        "missing_required_sites": missing_sites,
        "placeholder_row_count": len(placeholder_rows),
        "duplicate_row_count": len(duplicates),
        "schema_errors": schema_errors,
        "errors": errors,
        "warnings": warnings,
        "observed_site_count": len(observed_sites),
        "output_row_count": len(output_rows),
    }


def evaluate_one(expected, curated_rows, output_root):
    curated_row = curated_rows.get(normalize_key(expected.query))
    files = find_output_files(output_root, expected.output_dir)
    output_rows, schema_errors = collect_output_rows(files)

    identity = identity_result(expected, curated_row, output_rows)
    if identity["critical_identity_failure"]:
        ptm = {
            "ptm_score": 0,
            "missing_required_sites": list(expected.required_sites),
            "placeholder_row_count": 0,
            "duplicate_row_count": 0,
            "schema_errors": schema_errors,
            "errors": ["PTM scoring skipped because identity validation failed."],
            "warnings": [],
            "observed_site_count": 0,
            "output_row_count": len(output_rows),
        }
    else:
        ptm = ptm_result(expected, output_rows, schema_errors)

    final_score = identity["identity_score"] * ptm["ptm_score"] * 100
    errors = identity["errors"] + ptm["errors"]
    warnings = identity["warnings"] + ptm["warnings"]
    result = {
        "query": expected.query,
        "expected_protein_id": expected.expected_protein_id,
        "expected_organism": expected.expected_organism,
        "expected_uniprot": expected.expected_uniprot,
        "identity_score": identity["identity_score"],
        "ptm_score": ptm["ptm_score"],
        "final_score": final_score,
        "critical_identity_failure": identity["critical_identity_failure"],
        "missing_required_sites": ptm["missing_required_sites"],
        "placeholder_row_count": ptm["placeholder_row_count"],
        "duplicate_row_count": ptm["duplicate_row_count"],
        "schema_errors": ptm["schema_errors"],
        "output_files": [str(path) for path in files],
        "output_row_count": ptm["output_row_count"],
        "errors": errors,
        "warnings": warnings,
    }
    failure_class, signals, action = classify_result(result)
    result["failure_class"] = failure_class
    result["failure_signals"] = signals
    result["recommended_action"] = action
    return result


def classification_summary(results):
    counts = Counter(result["failure_class"] for result in results)
    total = len(results)
    valid_count = counts.get("VALID_SCRAPE", 0)
    identity_passes = sum(1 for result in results if result["identity_score"] == 1)
    ptm_complete = sum(
        1
        for result in results
        if result["identity_score"] == 1 and result["ptm_score"] == 1
    )
    blocked_like = (
        counts.get("CLOUDFLARE_LIKELY", 0)
        + counts.get("EMPTY_OUTPUT", 0)
        + counts.get("SCHEMA_FAILURE", 0)
    )
    return {
        "counts": dict(sorted(counts.items())),
        "usable_scrape_rate": valid_count / total if total else 0,
        "identity_pass_rate": identity_passes / total if total else 0,
        "ptm_completeness_rate": ptm_complete / identity_passes if identity_passes else 0,
        "blocked_or_invalid_rate": blocked_like / total if total else 0,
        "cloudflare_likely_rate": counts.get("CLOUDFLARE_LIKELY", 0) / total if total else 0,
    }


def evaluate(manifest_path, curated_ids_path, output_root):
    expected = load_expected_manifest(manifest_path)
    curated_rows = load_curated_ids(curated_ids_path)
    results = [evaluate_one(item, curated_rows, output_root) for item in expected]
    wrong_identity = sum(1 for result in results if result["critical_identity_failure"])
    wrong_protein_rate = wrong_identity / len(results) if results else 0
    final_score = sum(result["final_score"] for result in results) / len(results) if results else 0
    errors = [error for result in results for error in result["errors"]]
    warnings = [warning for result in results for warning in result["warnings"]]
    return {
        "status": "failed" if errors else "passed_with_warnings" if warnings else "passed",
        "protein_count": len(results),
        "wrong_protein_rate": wrong_protein_rate,
        "final_score": final_score,
        "classification": classification_summary(results),
        "errors": errors,
        "warnings": warnings,
        "results": results,
    }


def build_parser():
    parser = argparse.ArgumentParser(description="Identity-first deterministic evaluator for scraper outputs.")
    parser.add_argument("--manifest", required=True, help="Expected protein/PTM manifest CSV.")
    parser.add_argument("--curated-ids", required=True, help="Resolved protein ID CSV produced by the scraper.")
    parser.add_argument("--output-root", required=True, help="Root containing scraper output folders.")
    parser.add_argument("--json-output", help="Optional path for full JSON report.")
    parser.add_argument("--fail-on-error", action="store_true", help="Exit nonzero if evaluation status is failed.")
    return parser


def main():
    args = build_parser().parse_args()
    summary = evaluate(
        Path(args.manifest),
        Path(args.curated_ids),
        Path(args.output_root),
    )
    print(f"EVAL: status={summary['status']}", flush=True)
    print(f"EVAL: protein_count={summary['protein_count']}", flush=True)
    print(f"EVAL: wrong_protein_rate={summary['wrong_protein_rate']:.3f}", flush=True)
    print(f"EVAL: final_score={summary['final_score']:.1f}", flush=True)
    for error in summary["errors"][:20]:
        print(f"ERROR: {error}", flush=True)

    if args.json_output:
        output_path = Path(args.json_output)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w", encoding="utf-8") as handle:
            json.dump(summary, handle, indent=2, sort_keys=True)
        print(f"EVAL: wrote JSON report to {output_path}", flush=True)

    print("EVAL_JSON: " + json.dumps(summary, sort_keys=True), flush=True)
    if args.fail_on_error and summary["status"] == "failed":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
