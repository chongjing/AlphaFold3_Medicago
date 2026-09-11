#!/usr/bin/env python3
"""Resumable curl FTP upload of the 42 FASTQs to webin2.ebi.ac.uk.

Remote path (webin-cli convention): webin-cli/reads/<LIBRARY_NAME>/<FASTQ>
State file upload_state.tsv skips files already done with matching size.
Credentials: ENA_WEBIN_USER, ENA_WEBIN_PASSWORD.
"""
from __future__ import annotations

import csv
import os
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent
STATE = ROOT / "upload_state.tsv"
LOG = ROOT / "upload_fastq.log"
HOST = "ftp://webin2.ebi.ac.uk"
FIELDS = ["library_name", "filename", "local_bytes", "status", "utc", "note"]


def log(msg: str) -> None:
    line = f"{datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')} {msg}"
    print(line, flush=True)
    with LOG.open("a") as fh:
        fh.write(line + "\n")


def load_state() -> dict[str, dict]:
    if not STATE.exists():
        return {}
    with STATE.open() as fh:
        return {r["filename"]: r for r in csv.DictReader(fh, delimiter="\t")}


def write_state(state: dict) -> None:
    tmp = STATE.with_suffix(".tmp")
    with tmp.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS, delimiter="\t")
        w.writeheader()
        for name in sorted(state):
            w.writerow(state[name])
    tmp.replace(STATE)


def creds() -> tuple[str, str]:
    u = os.environ.get("ENA_WEBIN_USER")
    p = os.environ.get("ENA_WEBIN_PASSWORD")
    if not u or not p:
        raise SystemExit("set ENA_WEBIN_USER and ENA_WEBIN_PASSWORD")
    return u, p


def jobs() -> list[dict]:
    rows = []
    with (ROOT / "metadata.tsv").open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            for side in ("1", "2"):
                rows.append({
                    "library_name": row["library_name"],
                    "filename": row[f"fq{side}"],
                    "local_bytes": int(row[f"fq{side}_bytes"]),
                })
    return rows


def main() -> int:
    user, password = creds()
    state = load_state()
    todo = jobs()
    failed = 0
    for j in todo:
        lib, fn = j["library_name"], j["filename"]
        local = ROOT / fn
        if not local.exists():
            log(f"MISSING {fn}")
            failed += 1
            continue
        size = local.stat().st_size
        prev = state.get(fn)
        if prev and prev.get("status") == "done" and int(prev.get("local_bytes") or 0) == size:
            log(f"SKIP {lib} {fn}")
            continue
        remote = f"{HOST}/webin-cli/reads/{lib}/{fn}"
        ok = False
        for attempt in range(1, 4):
            t0 = time.time()
            log(f"PUT {lib} {fn} bytes={size} attempt={attempt}")
            cmd = ["curl", "-f", "-T", str(local),
                   "--retry", "5", "--retry-delay", "30",
                   "--connect-timeout", "60", "--max-time", "10800",
                   "--ftp-create-dirs", "--user", f"{user}:{password}", remote]
            try:
                subprocess.run(cmd, check=True)
                mbps = (size / 1e6) / (time.time() - t0) if time.time() - t0 else 0
                state[fn] = {"library_name": lib, "filename": fn,
                             "local_bytes": str(size), "status": "done",
                             "utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
                             "note": f"{mbps:.1f}MB/s"}
                write_state(state)
                ok = True
                break
            except subprocess.CalledProcessError as e:
                log(f"RETRY {lib} {fn} {e}")
                time.sleep(30)
        if not ok:
            state[fn] = {"library_name": lib, "filename": fn,
                         "local_bytes": str(size), "status": "failed",
                         "utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
                         "note": ""}
            write_state(state)
            failed += 1
            log(f"FAILED {lib} {fn}")
    done = sum(1 for j in todo if state.get(j["filename"], {}).get("status") == "done")
    log(f"finish done={done}/{len(todo)} failed={failed}")
    return 1 if failed or done != len(todo) else 0


if __name__ == "__main__":
    sys.exit(main())