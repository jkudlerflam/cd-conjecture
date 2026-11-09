"""C182: chemical potential mu(N) = E(N) - E(N-1) from local/eos (open chains, L = 40, phi = pi/2)."""
import json, glob, numpy as np, re
from pairband import valleys
E = {}
for f in glob.glob('local/eos/e_*.jsonl'):
    r = [json.loads(l) for l in open(f) if l.strip()]
    if not r or r[-1]['sweep'] != len(r[-1]['cfg']['chis'])-1: continue
    m = re.match(r'e_w([\d.]+)_N(\d+)_L40', r[0]['tag'])
    E[(float(m.group(1)), int(m.group(2)))] = (r[-1]['E'], r[-1]['trunc'])
for wa in sorted(set(k[0] for k in E)):
    em = valleys(1j*wa)[0][1]
    Ns = sorted(N for (w, N) in E if w == wa)
    print(f"|w|={wa}  eps_min={em:.4f}  (E(N)-N eps_min)/N and mu(N) = E(N)-E(N-1); mu2(N) = [E(N)-E(N-2)]/2 on even N")
    for N in Ns:
        e, tr = E[(wa, N)]
        s = f"  N={N:2d} n={N/40:.3f} e/N-eps={e/N-em:+.5f}"
        if (wa, N-1) in E: s += f" mu={e-E[(wa,N-1)][0]:+.5f}"
        if N % 2 == 0 and (wa, N-2) in E: s += f" mu2={(e-E[(wa,N-2)][0])/2:+.5f}"
        s += f" tr={tr:.0e}"
        print(s)
