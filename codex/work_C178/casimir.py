"""P4 (Casimir) and P3 (velocities) from twist-relaxed ground states.

For runs at the same (model, phi, |w|, n): E0_rel(L)/L = e - A / L^2 with A = pi (c v)_eff / 6, fitted from the
even-N sizes (two sizes: exact solve; three: least squares, residual printed). Velocities:
  v_J = pi L E''(theta) (flux stiffness of the relaxed ground-state branch),
  v_N = (L / pi) [E(N+2) + E(N-2) - 2 E(N)] / 4  (compressibility, relaxed energies; needs the N +- 2 runs),
  v_c = sqrt(v_J v_N)  (normalisation-free charge velocity for one or two components, see PROOF 4),
  v_k1 = L D_kappa1 / 2 pi (lowest M = 0 level at pseudospin momentum 2 pi / N, or the phonon for rods/legs).
Usage: python3 casimir.py files.jsonl"""
import sys, json, numpy as np
from collections import defaultdict
import analyze2 as A

recs = []
for fn in sys.argv[1:]:
    for line in open(fn):
        recs.append(json.loads(line))
summ = {}
for r in recs:
    o = A.summarize(r, verbose=False)
    summ[(r['model'], round(r['phi'], 4), round(r['wabs'], 4), r['L'], r['N'])] = o

groups = defaultdict(list)
for (model, phi, wa, L, N), o in summ.items():
    from fractions import Fraction
    groups[(model, phi, wa, Fraction(N, L))].append((L, N, o))

for key in sorted(groups, key=lambda k: (k[0], k[1], k[2], float(k[3]))):
    model, phi, wa, n = key
    rows = sorted(groups[key])
    print(f"\n{model} phi={phi} |w|={wa} n={n}:")
    for (L, N, o) in rows:
        vN = None
        lo = summ.get((model, phi, wa, L, N - 2)); hi = summ.get((model, phi, wa, L, N + 2))
        if lo and hi:
            vN = (L / np.pi) * (hi['E0'] + lo['E0'] - 2 * o['E0']) / 4
        vJ = o.get('vJ')
        vc = np.sqrt(vJ * vN) if (vN and vJ and vN > 0) else None
        if model in ('rods', 'legs'):
            k1 = o.get('rods_k1'); vk1 = L * k1[0] / (2 * np.pi) if k1 else None
        else:
            vk1 = L * o['D_k1'] / (2 * np.pi) if o.get('D_k1') else None
        o['vN'], o['vc'], o['vk1'] = vN, vc, vk1
        f = lambda x: f"{x:.3f}" if x is not None else "  -  "
        print(f"   L={L:3d} N={N:2d}: E0_rel={o['E0']:.6f} e_L={o['E0']/L:.6f}  vJ={f(vJ)} vN={f(vN)} v_c={f(vc)} v_k1={f(vk1)}" +
              (f" L*D_ps={L*o['D_ps']:.3f}" if o.get('D_ps') else ""))
    ev = [(L, o['E0']) for (L, N, o) in rows if N % 2 == 0]
    if len(ev) >= 2:
        Ls = np.array([x[0] for x in ev], float); y = np.array([x[1] for x in ev]) / Ls
        M = np.vstack([np.ones_like(Ls), -1 / Ls ** 2]).T
        sol, res, *_ = np.linalg.lstsq(M, y, rcond=None)
        cv = 6 * sol[1] / np.pi
        print(f"   Casimir (even N, L={list(map(int, Ls))}): e = {sol[0]:.6f}, (c v)_eff = {cv:.3f}" +
              (f", residual {np.sqrt(res[0]):.1e}" if len(res) else ""))
