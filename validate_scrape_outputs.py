import argparse
import csv
import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parent
CURATED_DIR = ROOT / "curated_protein_ids"
CURATED_ID_CSV = CURATED_DIR / "resolved_protein_ids.csv"
IGNORED_DIRS = {
    ".git",
    ".git_broken",
    ".venv",
    ".venv312",
    "__pycache__",
    "curated_protein_ids",
}


def is_ignored(path):
    return any(part in IGNORED_DIRS for part in path.parts)


def output_csv_files(root):
    files = []
    for path in root.rglob("*.csv"):
        relative = path.relative_to(root)
        if is_ignored(relative):
            continue
        files.append(path)
    return sorted(files, key=lambda item: str(item.relative_to(root)).lower())


def validate_curated_ids(path):
    result = {
        "path": str(path.relative_to(ROOT)) if path.exists() else str(path),
        "exists": path.exists(),
        "rows": 0,
        "errors": [],
        "warnings": [],
    }
    if not path.exists():
        result["warnings"].append("Curated protein ID cache does not exist yet.")
        return result

    try:
        with open(path, "r", encoding="utf-8", errors="replace", newline="") as handle:
            reader = csv.DictReader(handle)
            columns = set(reader.fieldnames or [])
            rows = list(reader)
    except Exception as exc:
        result["errors"].append(f"Could not read curated ID cache: {exc}")
        return result

    result["rows"] = len(rows)
    required = {"protein_name", "protein_id", "url", "organism", "source"}
    missing = sorted(required - columns)
    if missing:
        result["errors"].append(f"Missing curated ID columns: {', '.join(missing)}")
        return result

    for index, row in enumerate(rows):
        protein_name = str(row.get("protein_name", "")).strip()
        protein_id = row.get("protein_id")
        url = str(row.get("url", "")).strip()
        organism = str(row.get("organism", "")).strip().lower()

        if not protein_name:
            result["errors"].append(f"Row {index + 2}: missing protein_name")
        if not str(protein_id).split(".", 1)[0].isdigit():
            result["errors"].append(f"Row {index + 2}: invalid protein_id={protein_id}")
        if "proteinAction.action?id=" not in url:
            result["errors"].append(f"Row {index + 2}: URL is not a proteinAction URL")
        elif str(int(float(protein_id))) not in url:
            result["warnings"].append(f"Row {index + 2}: protein_id does not appear in URL")
        if organism and organism != "human":
            result["warnings"].append(f"Row {index + 2}: organism is '{organism}', expected human")

    return result


def validate_output_csv(path, root, deep=False):
    relative = path.relative_to(root)
    result = {
        "file": str(relative),
        "rows": None,
        "columns": [],
        "size_bytes": 0,
        "errors": [],
        "warnings": [],
    }
    try:
        result["size_bytes"] = path.stat().st_size
    except Exception as exc:
        result["errors"].append(f"Could not stat CSV: {exc}")
        return result

    if result["size_bytes"] == 0:
        result["warnings"].append("CSV file is empty.")
        return result

    if not deep:
        return result

    try:
        with open(path, "r", encoding="utf-8", errors="replace", newline="") as handle:
            reader = csv.DictReader(handle)
            columns = reader.fieldnames or []
            rows = list(reader)
    except Exception as exc:
        result["errors"].append(f"Could not read CSV: {exc}")
        return result

    result["rows"] = len(rows)
    result["columns"] = columns
    if not rows:
        result["warnings"].append("CSV has zero data rows.")
        return result

    if len(columns) < 2:
        result["warnings"].append("CSV has fewer than two columns.")

    parent_protein = path.parent.name.strip().lower()
    if "Protein" in columns:
        observed = {
            str(row.get("Protein", "")).strip().lower()
            for row in rows
            if str(row.get("Protein", "")).strip()
        }
        if observed and parent_protein not in observed:
            result["warnings"].append(
                f"Parent folder '{path.parent.name}' not found in Protein column values: {sorted(observed)[:5]}"
            )

    if "PubMed ID" in columns:
        pubmed_values = [str(row.get("PubMed ID", "")).strip() for row in rows if str(row.get("PubMed ID", "")).strip()]
        invalid_pubmed = [
            value
            for value in pubmed_values
            if value.lower() != "nan" and not re.fullmatch(r"\d{7,}", value.strip())
        ]
        if invalid_pubmed:
            result["warnings"].append(f"Found {len(invalid_pubmed)} non-numeric PubMed ID value(s).")

    return result


def build_summary(root, deep_csv=False, max_files=None, include_files=False):
    curated = validate_curated_ids(CURATED_ID_CSV)
    files = output_csv_files(root)
    if max_files is not None:
        files = files[:max_files]
    output_results = []
    for index, path in enumerate(files, start=1):
        if index == 1 or index % 100 == 0:
            print(f"VALIDATION: scanning output CSV {index}/{len(files)}", flush=True)
        output_results.append(validate_output_csv(path, root, deep=deep_csv))
    errors = curated["errors"] + [
        f"{item['file']}: {error}"
        for item in output_results
        for error in item["errors"]
    ]
    warnings = curated["warnings"] + [
        f"{item['file']}: {warning}"
        for item in output_results
        for warning in item["warnings"]
    ]
    total_rows = sum(item["rows"] for item in output_results if isinstance(item["rows"], int))
    proteins = sorted({path.parent.name for path in files})
    return {
        "status": "failed" if errors else "passed_with_warnings" if warnings else "passed",
        "curated_ids": curated,
        "output_csv_count": len(files),
        "output_row_count": total_rows,
        "deep_csv": deep_csv,
        "max_files": max_files,
        "protein_count": len(proteins),
        "proteins": proteins,
        "errors": errors,
        "warnings": warnings,
        "files": output_results if include_files else [],
    }


def build_parser():
    parser = argparse.ArgumentParser(description="Validate curated IDs and scraped PhosphoSitePlus CSV outputs.")
    parser.add_argument("--root", default=str(ROOT), help="Repository/output root to validate.")
    parser.add_argument("--json-output", help="Optional path to write full validation JSON.")
    parser.add_argument("--fail-on-warnings", action="store_true", help="Exit nonzero when warnings are present.")
    parser.add_argument("--deep-csv", action="store_true", help="Read every CSV and validate row-level content.")
    parser.add_argument("--max-files", type=int, help="Validate at most this many output CSV files.")
    parser.add_argument("--include-files", action="store_true", help="Include per-file validation details in JSON output.")
    return parser


def main():
    args = build_parser().parse_args()
    root = Path(args.root).resolve()
    summary = build_summary(
        root,
        deep_csv=args.deep_csv,
        max_files=args.max_files,
        include_files=args.include_files,
    )

    print(f"VALIDATION: status={summary['status']}", flush=True)
    print(f"VALIDATION: output_csv_count={summary['output_csv_count']}", flush=True)
    print(f"VALIDATION: output_row_count={summary['output_row_count']}", flush=True)
    print(f"VALIDATION: protein_count={summary['protein_count']}", flush=True)
    for error in summary["errors"][:20]:
        print(f"ERROR: {error}", flush=True)
    for warning in summary["warnings"][:20]:
        print(f"WARNING: {warning}", flush=True)

    if args.json_output:
        output_path = Path(args.json_output)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w", encoding="utf-8") as handle:
            json.dump(summary, handle, indent=2, sort_keys=True)
        print(f"VALIDATION: wrote JSON report to {output_path}", flush=True)

    print("VALIDATION_JSON: " + json.dumps(summary, sort_keys=True), flush=True)
    if summary["errors"] or (args.fail_on_warnings and summary["warnings"]):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
