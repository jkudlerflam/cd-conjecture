"""C194 (unregistered check): energies extrapolated linearly in the truncation error over the last two chi stages
(last sweep of each), then Delta_pb, then the registered L-fit on those.  Usage: python3 -E -s gap_xtrap.py [dirs]"""
import json, glob, os, sys
import numpy as np
from fitgap import fit_gap
st = {}
for d in sys.argv[1:] or ['out_main']:
    for f in glob.glob(os.path.join(d, '*.jsonl')):
        recs = [json.loads(l) for l in open(f) if l.strip()]
        fin = [r for r in recs if r.get('kind') == 'final']
        if not fin: continue
        c = fin[0]['cfg']; sw = [r for r in recs if r.get('kind') == 'sweep']
        stage = {}
        for r in sw: stage[r['chi']] = (r['E'], r['trunc'])
        chis = sorted(stage)[-2:]
        (E1, t1), (E2, t2) = stage[chis[0]], stage[chis[1]]
        Ex = E2 - t2 * (E1 - E2) / (t1 - t2) if t1 != t2 else E2
        k = (c['nom'], c['L'], c['N'])
        if k not in st or chis[1] > st[k][2]:
            st[k] = (E2, Ex, chis[1])
rows = {}
for (nom, L, N), (E, Ex, chi) in st.items():
    if N % 2: continue
    m, p = st.get((nom, L, N - 1)), st.get((nom, L, N + 1))
    if not (m and p): continue
    rows.setdefault(nom, []).append((L, p[0] + m[0] - 2 * E, p[1] + m[1] - 2 * Ex, chi))
for nom in sorted(rows, key=eval):
    r = sorted(rows[nom])
    print(f"n={nom}: " + '  '.join(f"L={L}(chi{c}) L*D={L*d:.4f} xtrap {L*dx:.4f}" for L, d, dx, c in r))
    if len(r) >= 3:
        f0 = fit_gap([x[0] for x in r], [x[1] for x in r]); f1 = fit_gap([x[0] for x in r], [x[2] for x in r])
        print(f"   D_inf raw {f0['D']:+.5f} +- {f0['err']:.5f};  truncation-extrapolated {f1['D']:+.5f} +- {f1['err']:.5f}")
