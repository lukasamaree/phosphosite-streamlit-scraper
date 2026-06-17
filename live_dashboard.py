import json
import subprocess
import sys
import tempfile
import threading
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs


ROOT = Path(__file__).resolve().parent
AGENTIC_WORKFLOW = ROOT / "agentic_phospho_workflow.py"
SCRAPER = ROOT / "phospho_group_scraper.py"
CURATION_DIR = ROOT / "curated_protein_ids"
LOOKUP_CSV = CURATION_DIR / "resolved_protein_ids.csv"
STATE_JSON = CURATION_DIR / "lookup_state.json"

CURATION_DIR.mkdir(exist_ok=True)

JOB = {
    "running": False,
    "returncode": None,
    "lines": [],
    "command": [],
}
LOCK = threading.Lock()


def parse_names(text):
    names = []
    normalized = text.replace("[", "\n").replace("]", "\n").replace(",", "\n")
    for line in normalized.splitlines():
        value = line.split("#", 1)[0].strip().strip("'\" \t")
        if value:
            names.append(value)
    return names


def append_log(line):
    with LOCK:
        JOB["lines"].append(line)
        JOB["lines"] = JOB["lines"][-500:]


def run_process(command):
    with LOCK:
        JOB["running"] = True
        JOB["returncode"] = None
        JOB["command"] = command
        JOB["lines"] = ["> " + " ".join(command)]

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
        clean = line.rstrip()
        if clean:
            append_log(clean)

    returncode = process.wait()
    append_log(f"Process exited with code {returncode}")
    with LOCK:
        JOB["running"] = False
        JOB["returncode"] = returncode


def start_agentic_lookup(proteins, delay, attempts):
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False, encoding="utf-8") as handle:
        for protein in proteins:
            handle.write(f"{protein}\n")
        names_file = handle.name

    command = [
        sys.executable,
        "-u",
        "-B",
        str(AGENTIC_WORKFLOW),
        "--protein-names-file",
        names_file,
        "--output-csv",
        str(LOOKUP_CSV),
        "--state-json",
        str(STATE_JSON),
        "--attempts",
        str(attempts),
        "--delay",
        str(delay),
    ]
    threading.Thread(target=run_process, args=(command,), daemon=True).start()


def start_scrape(delay):
    if not LOOKUP_CSV.exists():
        append_log("No resolved ID CSV found. Run agentic lookup first.")
        return

    ids = []
    with open(LOOKUP_CSV, "r", encoding="utf-8") as handle:
        header = handle.readline().strip().split(",")
        protein_id_index = header.index("protein_id")
        for line in handle:
            parts = line.strip().split(",")
            if len(parts) > protein_id_index and parts[protein_id_index].isdigit():
                ids.append(parts[protein_id_index])

    ids_file = CURATION_DIR / "protein_ids.txt"
    ids_file.write_text("\n".join(ids) + "\n", encoding="utf-8")
    command = [
        sys.executable,
        "-u",
        "-B",
        str(SCRAPER),
        "--protein-ids-file",
        str(ids_file),
        "--continue-on-error",
        "--delay",
        str(delay),
    ]
    threading.Thread(target=run_process, args=(command,), daemon=True).start()


PAGE = """<!doctype html>
<html>
<head>
  <meta charset="utf-8">
  <title>PhosphoSitePlus Agent Dashboard</title>
  <style>
    body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; margin: 0; background: #f7f8fa; color: #20242a; }
    main { max-width: 1180px; margin: 0 auto; padding: 28px; }
    h1 { margin: 0 0 6px; font-size: 28px; }
    .caption { color: #576170; margin-bottom: 22px; }
    .grid { display: grid; grid-template-columns: 420px 1fr; gap: 18px; align-items: start; }
    section { background: white; border: 1px solid #dde2e8; border-radius: 8px; padding: 18px; }
    textarea { width: 100%; min-height: 240px; box-sizing: border-box; font: 14px ui-monospace, SFMono-Regular, Menlo, monospace; }
    label { display: block; font-weight: 650; margin: 12px 0 6px; }
    input { width: 120px; padding: 7px; }
    button { padding: 9px 12px; border: 1px solid #315f9e; background: #386fb8; color: white; border-radius: 6px; cursor: pointer; margin: 12px 8px 0 0; }
    button.secondary { background: #ffffff; color: #2a4c78; }
    pre { background: #101720; color: #d9f2e6; padding: 14px; border-radius: 7px; min-height: 430px; max-height: 620px; overflow: auto; white-space: pre-wrap; }
    .status { font-weight: 700; margin-bottom: 10px; }
  </style>
</head>
<body>
<main>
  <h1>PhosphoSitePlus Agent Dashboard</h1>
  <div class="caption">Human-only workflow: search protein name, navigate to the human protein page, read the final proteinAction URL ID, checkpoint progress, retry failures.</div>
  <div class="grid">
    <section>
      <form id="runForm">
        <label>Protein symbols</label>
        <textarea name="proteins">TP53
AKT1
AKT2
AKT3
EGFR
MAPK1</textarea>
        <label>Delay between requests</label>
        <input name="delay" value="5" type="number" min="0" step="0.5">
        <label>Retry attempts</label>
        <input name="attempts" value="3" type="number" min="1" max="10">
        <br>
        <button name="action" value="lookup">Agentic Resolve IDs</button>
        <button class="secondary" name="action" value="scrape">Scrape Resolved IDs</button>
      </form>
    </section>
    <section>
      <div class="status" id="status">Idle</div>
      <pre id="log"></pre>
    </section>
  </div>
</main>
<script>
const form = document.getElementById('runForm');
form.addEventListener('submit', async (event) => {
  event.preventDefault();
  const submitter = event.submitter;
  const body = new URLSearchParams(new FormData(form));
  body.set('action', submitter.value);
  await fetch('/run', { method: 'POST', body });
});

async function refreshLog() {
  const response = await fetch('/logs');
  const data = await response.json();
  document.getElementById('status').textContent = data.running ? 'Running' : 'Idle';
  document.getElementById('log').textContent = data.lines.join('\\n');
}
setInterval(refreshLog, 1000);
refreshLog();
</script>
</body>
</html>
"""


class Handler(BaseHTTPRequestHandler):
    def do_GET(self):
        if self.path == "/logs":
            with LOCK:
                payload = json.dumps(JOB).encode("utf-8")
            self.send_response(200)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(payload)))
            self.end_headers()
            self.wfile.write(payload)
            return

        payload = PAGE.encode("utf-8")
        self.send_response(200)
        self.send_header("Content-Type", "text/html; charset=utf-8")
        self.send_header("Content-Length", str(len(payload)))
        self.end_headers()
        self.wfile.write(payload)

    def do_POST(self):
        length = int(self.headers.get("Content-Length", "0"))
        data = parse_qs(self.rfile.read(length).decode("utf-8"))
        action = data.get("action", ["lookup"])[0]
        delay = float(data.get("delay", ["5"])[0])
        attempts = int(data.get("attempts", ["3"])[0])

        with LOCK:
            running = JOB["running"]
        if running:
            self.send_response(409)
            self.end_headers()
            return

        if action == "scrape":
            start_scrape(delay)
        else:
            proteins = parse_names(data.get("proteins", [""])[0])
            start_agentic_lookup(proteins, delay, attempts)

        self.send_response(204)
        self.end_headers()

    def log_message(self, format, *args):
        return


if __name__ == "__main__":
    server = ThreadingHTTPServer(("127.0.0.1", 8512), Handler)
    print("Live dashboard running at http://127.0.0.1:8512")
    server.serve_forever()
