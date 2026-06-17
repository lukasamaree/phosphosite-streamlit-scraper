import argparse
import asyncio
import csv
import json
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

from phospho_group_scraper import parse_names_file, resolve_protein_name


ROOT = Path(__file__).resolve().parent
CURATION_DIR = ROOT / "curated_protein_ids"


@dataclass
class AgenticLookupConfig:
    attempts: int = 3
    delay: float = 5.0
    backoff: float = 2.0
    organism: str = "human"


def now_iso():
    return datetime.now(timezone.utc).isoformat()


def load_state(path):
    if not path.exists():
        return {"runs": {}, "updated_at": None}
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def save_state(path, state):
    path.parent.mkdir(parents=True, exist_ok=True)
    state["updated_at"] = now_iso()
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(state, handle, indent=2, sort_keys=True)


def write_lookup_csv(path, runs):
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    for protein_name, record in runs.items():
        if record.get("status") == "resolved":
            rows.append(
                {
                    "protein_name": protein_name,
                    "protein_id": record["protein_id"],
                    "url": record["url"],
                    "organism": record.get("organism", "human"),
                    "source": record.get("source", "agentic_lookup"),
                    "attempts": record.get("attempts", 1),
                    "status": record["status"],
                }
            )

    fieldnames = ["protein_name", "protein_id", "url", "organism", "source", "attempts", "status"]
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f"CHECKPOINT: wrote {len(rows)} resolved row(s) to {path}", flush=True)


async def run_agentic_lookup(protein_names, output_csv, state_json, config):
    state_path = Path(state_json)
    output_path = Path(output_csv)
    print(f"START: agentic human protein ID lookup for {len(protein_names)} protein(s)", flush=True)
    print(f"CONFIG: attempts={config.attempts}, delay={config.delay}s, backoff={config.backoff}x", flush=True)
    print(f"STATE: loading checkpoint from {state_path}", flush=True)
    state = load_state(state_path)
    runs = state.setdefault("runs", {})
    already_resolved = sum(1 for name in protein_names if runs.get(name, {}).get("status") == "resolved")
    print(f"STATE: {already_resolved}/{len(protein_names)} requested protein(s) already resolved", flush=True)

    for index, protein_name in enumerate(protein_names, start=1):
        existing = runs.get(protein_name)
        if existing and existing.get("status") == "resolved":
            print(
                f"[{index}/{len(protein_names)}] SKIP {protein_name}: already resolved as "
                f"{existing.get('protein_id')} from {existing.get('url')}",
                flush=True,
            )
            continue

        print(f"[{index}/{len(protein_names)}] START {protein_name}: resolving human URL ID", flush=True)
        last_error = None
        for attempt in range(1, config.attempts + 1):
            try:
                print(
                    f"[{index}/{len(protein_names)}] {protein_name}: attempt {attempt}/{config.attempts} "
                    "launching browser search",
                    flush=True,
                )
                result = await resolve_protein_name(protein_name, config.organism)
                runs[protein_name] = {
                    **result,
                    "status": "resolved",
                    "attempts": attempt,
                    "resolved_at": now_iso(),
                }
                print(
                    f"[{index}/{len(protein_names)}] RESOLVED {protein_name}: "
                    f"id={result['protein_id']} url={result['url']}",
                    flush=True,
                )
                print(f"CHECKPOINT: saving state to {state_path}", flush=True)
                save_state(state_path, state)
                write_lookup_csv(output_path, runs)
                break
            except Exception as exc:
                last_error = str(exc)
                runs[protein_name] = {
                    "status": "retrying" if attempt < config.attempts else "failed",
                    "attempts": attempt,
                    "last_error": last_error,
                    "updated_at": now_iso(),
                }
                print(f"CHECKPOINT: saving failed attempt state to {state_path}", flush=True)
                save_state(state_path, state)
                write_lookup_csv(output_path, runs)
                print(f"ERROR: {protein_name} attempt {attempt}/{config.attempts} failed: {last_error}", flush=True)
                if attempt < config.attempts:
                    wait_seconds = config.delay * (config.backoff ** (attempt - 1))
                    print(f"WAIT: {protein_name} retrying after {wait_seconds:.1f}s", flush=True)
                    await asyncio.sleep(wait_seconds)

        if runs[protein_name].get("status") != "resolved":
            print(f"FAILED: {protein_name}: {last_error}", flush=True)

        if index < len(protein_names) and config.delay > 0:
            print(f"WAIT: sleeping {config.delay:.1f}s before next protein", flush=True)
            await asyncio.sleep(config.delay)

    write_lookup_csv(output_path, runs)
    resolved = sum(1 for record in runs.values() if record.get("status") == "resolved")
    failed = sum(1 for record in runs.values() if record.get("status") == "failed")
    print(f"DONE: agentic lookup complete: {resolved} resolved, {failed} failed", flush=True)
    return state


def build_parser():
    parser = argparse.ArgumentParser(description="Agentic PhosphoSitePlus human protein ID lookup.")
    parser.add_argument("--protein-names", nargs="+")
    parser.add_argument("--protein-names-file")
    parser.add_argument("--output-csv", default=str(CURATION_DIR / "resolved_protein_ids.csv"))
    parser.add_argument("--state-json", default=str(CURATION_DIR / "lookup_state.json"))
    parser.add_argument("--attempts", type=int, default=3)
    parser.add_argument("--delay", type=float, default=5.0)
    parser.add_argument("--backoff", type=float, default=2.0)
    return parser


def collect_names(args):
    names = []
    if args.protein_names:
        names.extend(args.protein_names)
    if args.protein_names_file:
        names.extend(parse_names_file(args.protein_names_file))
    return names


if __name__ == "__main__":
    args = build_parser().parse_args()
    protein_names = collect_names(args)
    if not protein_names:
        raise SystemExit("No protein names provided.")

    config = AgenticLookupConfig(
        attempts=args.attempts,
        delay=args.delay,
        backoff=args.backoff,
    )
    asyncio.run(run_agentic_lookup(protein_names, args.output_csv, args.state_json, config))
