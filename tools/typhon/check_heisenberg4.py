"""Check a fetched four-site Heisenberg spectrum against the exact singlet and triplet."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import tarfile


def check(run):
    job = (run / "job-id").read_text().strip().split(";")[0]
    runtime = json.loads((run / "logs/runtime.json").read_text())
    with (run / "accounting.psv").open() as f:
        rows = [r for r in csv.DictReader(f, delimiter="|") if r["JobIDRaw"] == job]
    if len(rows) != 1 or (rows[0]["State"], rows[0]["ExitCode"]) != ("COMPLETED", "0:0"):
        raise ValueError(f"Expected a completed allocation and exit 0:0 for job {job}")
    node = runtime["node"].split(".")[0]
    if runtime["job_id"] != job or rows[0]["NodeList"] != node:
        raise ValueError("Compute-node log and Slurm accounting identify different runs")

    record = json.loads((run / "results/spectrum.jsonl").read_text())
    if (record["L"] != 4 or not record["periodic"] or record["distance"] != 1
            or record["couplings"] != {"Jx": 0.25, "Jy": 0.25, "Jz": 0.25}
            or record["requested_sectors"] != [0, 1]):
        raise ValueError("Expected the periodic four-site Heisenberg ring in both parity sectors")
    with tarfile.open(run / "snapshot.tar.gz") as archive:
        for key, path in (("source_sha256", "contrib/d2_order_one_20261110/ed_screen.py"),
                          ("config_sha256", "tools/typhon/heisenberg4.json")):
            with archive.extractfile("src/" + path) as f:
                if hashlib.sha256(f.read()).hexdigest() != record["provenance"][key]:
                    raise ValueError(f"The solver used different bytes for {path}")

    # Both exact observables are O(1). Allow ten times the requested eigensolver tolerance.
    tol = 1e-11
    for name, exact in (("ground_energy", -2.0), ("ground_gap", 1.0)):
        value = record[name]
        if not math.isclose(value, exact, abs_tol=tol, rel_tol=0):
            raise ValueError(f"{name} = {value}, expected {exact} within {tol:g}")
        print(f"{name} = {value} (exact {exact:g}, absolute error {abs(value - exact):.3g})")
    residual = max(r for sector in record["sectors"] for r in sector["residuals"])
    if not residual <= tol:
        raise ValueError(f"Eigenpair residual {residual} exceeds {tol:g}")
    print(f"Job {job} on {node}: COMPLETED, exit 0:0, max residual {residual:.3g}, tolerance {tol:g}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path, help="downloaded run directory")
    check(parser.parse_args().run)
