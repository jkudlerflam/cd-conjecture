"""Check the determinant kernel against the sparse-resolvent kernel of C94 (inherited_second_order)."""
import sys, json, time
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / 'codex' / 'work_C94'))
from inherited_second_order import susceptibility
from kernel import kernel
rows = []
for N, V in [(4, 4), (5, 3), (6, 6), (7, 5), (8, 8), (9, 4), (30, 1), (10, 6)]:
    t0 = time.time()
    chi_ed, res = susceptibility(N, V)
    chi, d = kernel(N, V)
    err = float(np.max(np.abs(chi - chi_ed[0])))
    rows.append(dict(N=N, V=V, max_abs_diff=err, chi0=float(chi[0]), sum_rule_error=d['sum_rule_error'],
                     max_imag=d['max_imag'], ed_residual=res))
    print(rows[-1], round(time.time() - t0, 1), flush=True)
(HERE / 'check_kernel_ed.json').write_text(json.dumps(rows, indent=1) + '\n')
