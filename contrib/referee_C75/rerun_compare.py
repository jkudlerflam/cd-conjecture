"""(4) Compare outputs regenerated in a scratch copy of codex/work_C75 (after
running run_certificate.py, checks.py and selector_checks.py there) with the
committed files. Usage: python3 rerun_compare.py SCRATCH_WORK_C75 REPO_WORK_C75"""
import sys, json
import numpy as np
S, O = sys.argv[1], sys.argv[2]
out = {}
def load(p):
    d = json.load(open(p))
    if isinstance(d, dict): d.pop('seconds', None)
    return d
def close(a, b, tol=1e-9):
    if isinstance(a, dict): return a.keys() == b.keys() and all(close(a[k], b[k], tol) for k in a)
    if isinstance(a, list): return len(a) == len(b) and all(close(x, y, tol) for x, y in zip(a, b))
    if isinstance(a, float): return abs(a - b) <= tol * max(1, abs(a))
    return a == b
for f in ['certificate_n15.json', 'root_certificate_exact_n15.json', 'root_certificate_n15.json',
          'full_loop_bands_n15.json', 'rational_bands_n15.json', 'finite_budget.json',
          'finite_budget_exact.json', 'selector_certificate.json', 'small_checks.json']:
    a, b = load(f'{S}/{f}'), load(f'{O}/{f}')
    out[f] = 'identical' if a == b else ('equal within 1e-9 relative (float noise)' if close(a, b) else 'DIFFERENT')
for f in ['certificate_n15_D1000.npz', 'derivative_n15.npz', 'full_loop_n15.npz',
          'full_loop_whitened_n15.npz', 'whitened_n15.npz']:
    a, b = np.load(f'{S}/{f}', allow_pickle=True), np.load(f'{O}/{f}', allow_pickle=True)
    out[f] = 'identical' if set(a.files) == set(b.files) and all(np.array_equal(a[k], b[k]) for k in a.files) else 'DIFFERENT'
print(json.dumps(out, indent=1)); json.dump(out, open('rerun_compare.json', 'w'), indent=1)
