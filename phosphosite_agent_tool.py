import argparse
import asyncio
import csv
import json
from pathlib import Path

from agentic_phospho_workflow import AgenticLookupConfig, load_state, run_agentic_lookup
from phospho_group_scraper import parse_ids_file, parse_names_file, run_protein_batch
from scraper_eval.evaluate_scraper_outputs import evaluate


ROOT = Path(__file__).resolve().parent
CURATION_DIR = ROOT / "curated_protein_ids"
DEFAULT_STATE_JSON = CURATION_DIR / "lookup_state.json"
DEFAULT_LOOKUP_CSV = CURATION_DIR / "resolved_protein_ids.csv"
DEFAULT_IDS_TXT = CURATION_DIR / "protein_ids.txt"
DEFAULT_SCRAPE_STATE_JSON = CURATION_DIR / "scrape_state.json"


def unique_names(names):
    seen = set()
    unique = []
    for name in names:
        value = str(name).strip()
        key = value.upper()
        if value and key not in seen:
            seen.add(key)
            unique.append(value)
    return unique


def collect_names(args):
    names = []
    if getattr(args, "protein_names", None):
        names.extend(args.protein_names)
    if getattr(args, "protein_names_file", None):
        names.extend(parse_names_file(args.protein_names_file))
    return unique_names(names)


def collect_ids(args):
    ids = []
    if getattr(args, "protein_ids", None):
        ids.extend(args.protein_ids)
    if getattr(args, "protein_ids_file", None):
        ids.extend(parse_ids_file(args.protein_ids_file, "protein ID"))
    return list(dict.fromkeys(ids))


def read_lookup_csv(path):
    path = Path(path)
    if not path.exists():
        return []
    with open(path, "r", encoding="utf-8", errors="replace", newline="") as handle:
        return list(csv.DictReader(handle))


def row_protein_key(row):
    return str(row.get("protein_name", "")).strip().upper()


def row_protein_id(row):
    value = str(row.get("protein_id", "")).strip()
    if not value:
        return None
    return int(float(value))


def resolved_rows_for_names(names, lookup_csv):
    lookup_rows = read_lookup_csv(lookup_csv)
    if not lookup_rows or not names:
        return [], names

    requested = {name.upper() for name in names}
    found = [row for row in lookup_rows if row_protein_key(row) in requested]
    found_keys = {row_protein_key(row) for row in found}
    missing = [name for name in names if name.upper() not in found_keys]
    return found, missing


def write_ids_file(ids, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        for protein_id in ids:
            handle.write(f"{int(protein_id)}\n")
    return path


def maybe_evaluate_outputs(args):
    manifest = getattr(args, "eval_manifest", None)
    if not manifest:
        return None

    output_root = Path(getattr(args, "output_root", ROOT)).resolve()
    summary = evaluate(Path(manifest), Path(args.lookup_csv), output_root)
    report_path = getattr(args, "eval_report", None)
    if report_path:
        report_path = Path(report_path)
        report_path.parent.mkdir(parents=True, exist_ok=True)
        with open(report_path, "w", encoding="utf-8") as handle:
            json.dump(summary, handle, indent=2, sort_keys=True)
        print(f"EVAL: wrote scraper evaluation report to {report_path.resolve()}", flush=True)

    print(
        "EVAL: "
        f"status={summary['status']} "
        f"final_score={summary['final_score']:.1f} "
        f"wrong_protein_rate={summary['wrong_protein_rate']:.3f} "
        f"classification={summary['classification']['counts']}",
        flush=True,
    )
    return summary


def summarize_lookup_state(state_json, lookup_csv, requested_names=None):
    state = load_state(Path(state_json))
    runs = state.get("runs", {})
    resolved_records = {
        name: record
        for name, record in runs.items()
        if record.get("status") == "resolved"
    }
    failed_records = {
        name: record
        for name, record in runs.items()
        if record.get("status") == "failed"
    }
    saved_rows = read_lookup_csv(lookup_csv)
    summary = {
        "state_json": str(Path(state_json).resolve()),
        "lookup_csv": str(Path(lookup_csv).resolve()),
        "state_updated_at": state.get("updated_at"),
        "requested_names": requested_names or [],
        "resolved_count": len(resolved_records),
        "failed_count": len(failed_records),
        "csv_rows": len(saved_rows),
        "resolved": [
            {
                "protein_name": name,
                "protein_id": record.get("protein_id"),
                "url": record.get("url"),
                "organism": record.get("organism", "human"),
                "status": record.get("status"),
            }
            for name, record in resolved_records.items()
        ],
        "failed": [
            {
                "protein_name": name,
                "error": record.get("last_error"),
                "attempts": record.get("attempts"),
                "status": record.get("status"),
            }
            for name, record in failed_records.items()
        ],
    }
    return summary


async def resolve_ids(args):
    names = collect_names(args)
    if not names:
        raise SystemExit("No protein names provided.")

    saved_rows, missing_names = resolved_rows_for_names(names, args.lookup_csv)
    if args.only_missing and not missing_names:
        print(f"SKIP: all {len(names)} protein name(s) already have saved IDs", flush=True)
    else:
        names_to_resolve = missing_names if args.only_missing else names
        print(
            f"TOOL: resolving {len(names_to_resolve)} protein name(s); "
            f"{len(saved_rows)} saved ID(s) reused",
            flush=True,
        )
        config = AgenticLookupConfig(
            attempts=args.attempts,
            delay=args.delay,
            backoff=args.backoff,
            cloudflare_cooldown=args.cloudflare_cooldown,
            delay_jitter=args.delay_jitter,
            organism="human",
        )
        await run_agentic_lookup(names_to_resolve, args.lookup_csv, args.state_json, config)

    summary = summarize_lookup_state(args.state_json, args.lookup_csv, names)
    print("SUMMARY_JSON: " + json.dumps(summary, sort_keys=True), flush=True)
    return summary


async def scrape_ids(args):
    ids = collect_ids(args)
    names = collect_names(args)

    if names:
        saved_rows, missing_names = resolved_rows_for_names(names, args.lookup_csv)
        if missing_names and args.resolve_missing:
            print(f"TOOL: resolving {len(missing_names)} missing ID(s) before scrape", flush=True)
            config = AgenticLookupConfig(
                attempts=args.attempts,
                delay=args.delay,
                backoff=args.backoff,
                cloudflare_cooldown=args.cloudflare_cooldown,
                delay_jitter=args.delay_jitter,
                organism="human",
            )
            await run_agentic_lookup(missing_names, args.lookup_csv, args.state_json, config)
            saved_rows, missing_names = resolved_rows_for_names(names, args.lookup_csv)

        if missing_names:
            print(f"WARNING: missing saved IDs for: {', '.join(missing_names)}", flush=True)
        ids.extend(
            protein_id
            for protein_id in (row_protein_id(row) for row in saved_rows)
            if protein_id is not None
        )

    ids = list(dict.fromkeys(int(protein_id) for protein_id in ids))
    if not ids:
        raise SystemExit("No protein IDs available to scrape.")

    ids_file = write_ids_file(ids, args.ids_file)
    print(f"TOOL: scraping {len(ids)} protein ID(s) from {ids_file}", flush=True)
    scrape_result = await run_protein_batch(
        ids,
        args.delay,
        args.continue_on_error,
        args.delay_jitter,
        args.scrape_state,
        args.force_rescrape,
    )

    summary = {
        "protein_ids": ids,
        "ids_file": str(ids_file.resolve()),
        "lookup_csv": str(Path(args.lookup_csv).resolve()),
        "state_json": str(Path(args.state_json).resolve()),
        "scrape_status": scrape_result.get("status", "complete") if scrape_result else "complete",
        "failed_protein_ids": scrape_result.get("failed_protein_ids", []) if scrape_result else [],
        "failed_site_ids_by_protein": scrape_result.get("failed_site_ids_by_protein", {}) if scrape_result else {},
    }
    evaluation = maybe_evaluate_outputs(args)
    if evaluation is not None:
        summary["evaluation"] = evaluation
    print("SUMMARY_JSON: " + json.dumps(summary, sort_keys=True), flush=True)
    return summary


async def run_full_workflow(args):
    args.only_missing = True
    args.resolve_missing = True
    await resolve_ids(args)
    return await scrape_ids(args)


def add_shared_lookup_args(parser):
    parser.add_argument("--protein-names", nargs="+", help="Protein symbols/names to resolve.")
    parser.add_argument("--protein-names-file", help="Text file with one protein symbol/name per line.")
    parser.add_argument("--lookup-csv", default=str(DEFAULT_LOOKUP_CSV), help="Resolved ID CSV cache.")
    parser.add_argument("--state-json", default=str(DEFAULT_STATE_JSON), help="Resolved ID JSON checkpoint.")
    parser.add_argument("--attempts", type=int, default=3, help="Lookup attempts per missing protein.")
    parser.add_argument("--delay", type=float, default=5.0, help="Seconds between requests.")
    parser.add_argument("--delay-jitter", type=float, default=0.5, help="Randomize scraper waits by this fraction around --delay.")
    parser.add_argument("--backoff", type=float, default=2.0, help="Retry backoff multiplier.")
    parser.add_argument("--cloudflare-cooldown", type=float, default=60.0, help="Seconds to wait after Cloudflare challenges.")
    parser.add_argument("--scrape-state", default=str(DEFAULT_SCRAPE_STATE_JSON), help="Resolved-ID scrape checkpoint JSON.")
    parser.add_argument("--force-rescrape", action="store_true", help="Ignore completed site IDs in the scrape checkpoint.")


def add_evaluation_args(parser):
    parser.add_argument(
        "--eval-manifest",
        help="Optional expected manifest CSV. When supplied, evaluate and classify outputs after scraping.",
    )
    parser.add_argument(
        "--eval-report",
        help="Optional JSON path for the post-scrape evaluation/classification report.",
    )
    parser.add_argument(
        "--output-root",
        default=str(ROOT),
        help="Root containing scraper output folders. Default: repository root/current scraper output root.",
    )


def build_parser():
    parser = argparse.ArgumentParser(
        description="Reusable PhosphoSitePlus tool for agents: resolve human protein IDs and scrape from IDs."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    resolve_parser = subparsers.add_parser("resolve", help="Resolve and save human proteinAction IDs.")
    add_shared_lookup_args(resolve_parser)
    resolve_parser.add_argument(
        "--all",
        action="store_false",
        dest="only_missing",
        help="Resolve all requested proteins again instead of only missing ones.",
    )
    resolve_parser.set_defaults(only_missing=True)

    scrape_parser = subparsers.add_parser("scrape", help="Scrape protein data from saved or supplied IDs.")
    add_shared_lookup_args(scrape_parser)
    scrape_parser.add_argument("--protein-ids", nargs="+", type=int, help="ProteinAction IDs to scrape.")
    scrape_parser.add_argument("--protein-ids-file", help="Text file with one proteinAction ID per line.")
    scrape_parser.add_argument("--ids-file", default=str(DEFAULT_IDS_TXT), help="Scratch ID file written for scraping.")
    scrape_parser.add_argument(
        "--resolve-missing",
        action="store_true",
        help="When protein names are supplied, resolve missing saved IDs before scraping.",
    )
    scrape_parser.add_argument("--continue-on-error", action="store_true", help="Continue scraping after failures.")
    add_evaluation_args(scrape_parser)

    run_parser = subparsers.add_parser("run", help="Resolve missing IDs, save them, then scrape them.")
    add_shared_lookup_args(run_parser)
    run_parser.add_argument("--protein-ids", nargs="+", type=int, help="Extra proteinAction IDs to scrape.")
    run_parser.add_argument("--protein-ids-file", help="Extra proteinAction IDs file.")
    run_parser.add_argument("--ids-file", default=str(DEFAULT_IDS_TXT), help="Scratch ID file written for scraping.")
    run_parser.add_argument("--continue-on-error", action="store_true", help="Continue scraping after failures.")
    add_evaluation_args(run_parser)

    return parser


async def main():
    args = build_parser().parse_args()
    if args.command == "resolve":
        await resolve_ids(args)
    elif args.command == "scrape":
        await scrape_ids(args)
    elif args.command == "run":
        await run_full_workflow(args)


if __name__ == "__main__":
    asyncio.run(main())
