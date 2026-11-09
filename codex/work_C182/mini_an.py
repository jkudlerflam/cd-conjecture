"""C182: two-size summary of local DMRG runs (any directory of run.py jsonl outputs with _L<L> tags)."""
import numpy as np, json, glob, sys, os
from analyze import Sbar, c_prof
d = sys.argv[1]
runs = {}
for f in glob.glob(f'{d}/*.jsonl'):
    r = [json.loads(l) for l in open(f) if l.strip()]
    if r: runs[r[0]['tag']] = r
sers = {}
for t, r in runs.items():
    s, L = t.rsplit('_L', 1)
    if '_' in L: continue
    sers.setdefault(s, {})[int(L)] = r
for s in sorted(sers):
    Ls = sorted(sers[s]); last = {L: sers[s][L][-1] for L in Ls}
    out = f"{s:28s}"
    for L in Ls:
        r = last[L]; m = r.get('meas')
        out += f" | L={L} sw={r['sweep']+1}/{len(r['cfg']['chis'])} tr={r['trunc']:.0e} S1={Sbar(r['S'],L,1):.3f}"
        if m:
            nm = np.mean(m['dens'][L//4:3*L//4]); W = sorted([m['Wa_f'], m['Wb_f']])
            out += f" nmid={nm:.3f} Wr={W[0]/max(W[1],1e-12):.2f}"
    if len(Ls) >= 2:
        L1, L2 = Ls[0], Ls[-1]
        for W in (1, L1//8):
            c = 6*(Sbar(last[L2]['S'], L2, W)-Sbar(last[L1]['S'], L1, W))/np.log(L2/L1)
            out += f" | c({L1}->{L2},W={W})={c:.2f}"
        out += f" cprof({L2})={c_prof(last[L2]['S'], L2, sm=8):.2f}"
    print(out)
