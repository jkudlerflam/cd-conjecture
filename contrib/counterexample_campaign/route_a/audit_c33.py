"""Rebuild and independently rerun the C33 n=23 certificate in a temp directory.

Original collaborators' files are only read. Requires numpy and scipy.
Writes a small audit summary here; deletes its temporary copied workspace.
"""
from pathlib import Path
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SOURCE = ROOT / "codex/work_C33"


def run():
    with tempfile.TemporaryDirectory(prefix="cd-c33-independent-audit-") as td:
        temp = Path(td)
        for p in SOURCE.glob("*.py"):
            shutil.copy2(p, temp/p.name)
        # The derivative matrix is only a conditioning proposal. integer_ldl
        # checks its exact congruence identity; final tests use the new loop.
        for name in ("derivative_n23.npz", "certificate_n23.json"):
            shutil.copy2(SOURCE/name, temp/name)
        env = dict(os.environ, OPENBLAS_NUM_THREADS="1")
        logs = []
        for script in ("full_loop_integer.py", "certify_roots.py"):
            result = subprocess.run([sys.executable, script], cwd=temp,
                env=env, check=True, text=True, capture_output=True, timeout=210)
            logs.append(dict(script=script, stdout=result.stdout.strip()))
        x = np.load(SOURCE/"full_loop_n23.npz", allow_pickle=True)
        y = np.load(temp/"full_loop_n23.npz", allow_pickle=True)
        matches = {k: bool(np.array_equal(x[k], y[k])) for k in x.files}
        assert all(matches.values())
        cert = json.loads((temp/"root_certificate_n23.json").read_text())
        # Exact rational decisions are made in certify_roots.py. These are
        # outward rounded summaries, not replacement acceptance criteria.
        assert cert["root_count"] == 46
        assert cert["finite_g_power"] == -365
        result = dict(status="independently_rerun_and_mathematically_audited",
            generated_arrays_equal=matches, certificate=cert,
            script_sha256={name:hashlib.sha256((SOURCE/name).read_bytes()).hexdigest()
                for name in ("full_loop_integer.py", "certify_roots.py",
                             "rational_bands.py")},
            logs=logs)
        (HERE/"c33_audit.json").write_text(json.dumps(result, indent=2)+"\n")
        print(json.dumps({k:v for k,v in result.items() if k != "logs"}, indent=2))


if __name__ == "__main__":
    run()
