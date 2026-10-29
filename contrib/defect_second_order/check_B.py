"""Check B against C89's operator-valued second-order coefficients (sparse resolvent, full clock space)."""
import json
from pathlib import Path
from kernel import kernel
from sources import selector, B_coeff
HERE = Path(__file__).resolve().parent
ref = json.loads((HERE.parents[1] / 'codex' / 'work_C89' / 'second_order_results.json').read_text())
cache = {}; rows = []; seen = set()
for r in ref['rows']:
    key = (r['N'], r['V'], tuple(r['lengths']), tuple(r['pins']))
    if key in seen: continue
    seen.add(key)
    if (r['N'], r['V']) not in cache: cache[r['N'], r['V']] = kernel(r['N'], r['V'])[0]
    chi = cache[r['N'], r['V']]
    mu = selector(r['b'])
    good = [r['b']] * (r['N'] // r['b'])
    Bg = B_coeff(chi, good, [False] * len(good), mu)
    Bd = B_coeff(chi, r['lengths'], r['pins'], mu)
    rows.append(dict(N=r['N'], V=r['V'], lengths=r['lengths'], pins=r['pins'], B_good=Bg, B_defect=Bd,
                     C89_B_good=r['B_reference'], C89_B_defect=r['B_defect'],
                     diff=max(abs(Bg - r['B_reference']), abs(Bd - r['B_defect']))))
    print(rows[-1])
(HERE / 'check_B.json').write_text(json.dumps(rows, indent=1) + '\n')
