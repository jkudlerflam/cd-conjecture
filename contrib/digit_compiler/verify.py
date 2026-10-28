#!/usr/bin/env python3
"""Rerun every new component and integration check, preserving a manifest."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
CHECKS = (
    "alu/checks.py", "alu/bounds_checks.py", "memory/checks.py",
    "memory/vm_checks.py", "audit/checks.py", "audit/single_head_checks.py",
    "audit/reversible_head_checks.py", "audit/gate_checks.py",
    "checks.py", "generator_checks.py", "reversible_checks.py",
)


def main():
    source_hashes = {str(p.relative_to(HERE)): hashlib.sha256(p.read_bytes()).hexdigest()
                     for p in sorted(HERE.rglob("*.py"))}
    rows = []
    started = time.perf_counter()
    for name in CHECKS:
        before = time.perf_counter()
        proc = subprocess.run([sys.executable, str(HERE/name)], cwd=REPO,
                              capture_output=True, text=True)
        if proc.returncode:
            print(proc.stdout)
            print(proc.stderr, file=sys.stderr)
            raise RuntimeError(f"Check failed: {name}")
        row = dict(check=name, status="PASS", seconds=round(time.perf_counter()-before, 3),
                   stdout_sha256=hashlib.sha256(proc.stdout.encode()).hexdigest())
        rows.append(row)
        print(json.dumps(row), flush=True)
    assert source_hashes == {str(p.relative_to(HERE)): hashlib.sha256(p.read_bytes()).hexdigest()
                             for p in sorted(HERE.rglob("*.py"))}, "Source changed during verification"
    result = dict(status="PASS", checks=rows, source_sha256=source_hashes,
                  runtime_seconds=round(time.perf_counter()-started, 3))
    (HERE/"verification_summary.json").write_text(json.dumps(result, indent=2)+"\n")


if __name__ == "__main__":
    main()
