"""Rerun the C33 n=23 finite-loop generator and root certificate in a scratch copy.

Reads codex/work_C33 only; writes c33_rerun.json here.  This is a
reproducibility rerun of the collaborators' own code, not an independent
re-derivation of the certificate.
"""
from pathlib import Path
import json, os, shutil, subprocess, sys, tempfile, time
import numpy as np

HERE = Path(__file__).resolve().parent
SRC = HERE.parents[1] / 'codex' / 'work_C33'
with tempfile.TemporaryDirectory(prefix='referee-d24-c33-') as td:
    tmp = Path(td)
    for p in SRC.glob('*.py'):
        shutil.copy2(p, tmp / p.name)
    for name in ('derivative_n23.npz', 'certificate_n23.json'):
        shutil.copy2(SRC / name, tmp / name)
    env = dict(os.environ, OPENBLAS_NUM_THREADS='1')
    logs = []
    for script in ('full_loop_integer.py', 'certify_roots.py'):
        t0 = time.time()
        r = subprocess.run([sys.executable, '-B', script], cwd=tmp, env=env, check=True,
                           text=True, capture_output=True, timeout=1200)
        logs.append(dict(script=script, seconds=round(time.time() - t0, 1),
                         stdout_tail=r.stdout.strip()[-1500:]))
    a = np.load(SRC / 'full_loop_n23.npz', allow_pickle=True)
    b = np.load(tmp / 'full_loop_n23.npz', allow_pickle=True)
    same = {k: bool(np.array_equal(a[k], b[k])) for k in a.files}
    cert = json.loads((tmp / 'root_certificate_n23.json').read_text())
    keys = [k for k in cert if not isinstance(cert[k], (list, dict))]
    out = dict(arrays_identical_to_stored=same, certificate_scalars={k: cert[k] for k in keys},
               logs=logs)
    (HERE / 'c33_rerun.json').write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps(dict(arrays=same, cert={k: cert[k] for k in keys}), indent=1)[:3000])
