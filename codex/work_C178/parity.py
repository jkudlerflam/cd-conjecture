"""Odd-even (parity) gap and flux-period diagnostics from twist-relaxed ring energies.

Delta_pb(L, N) = E(N+1) + E(N-1) - 2 E(N) at even N, all energies twist-relaxed (min over flux and momentum).
Gapless two-component liquid: L * Delta_pb tends to a constant (charge and pseudospin zero modes, both ~ 1/L).
Pair liquid (pseudospin gap Delta_s): Delta_pb -> 2 Delta_s, so L * Delta_pb grows linearly with L.
Also printed: delta_pi = min_K E(theta = pi) - min_K E(theta = 0) at fixed N (flux-period test: a liquid of
charge-2 pairs has period pi in theta, so delta_pi -> 0 faster than 1/L; a two-component liquid has L*|delta_pi| ~ const).
Usage: python3 parity.py files.jsonl"""
import sys, json, numpy as np
from collections import defaultdict
import analyze2 as A

E = {}
dpi = {}
for fn in sys.argv[1:]:
    for line in open(fn):
        r = json.loads(line)
        key = (r['model'], round(r['phi'], 4), round(r['wabs'], 4), r['L'], r['N'])
        o = A.summarize(r, verbose=False)
        E[key] = min(E.get(key, 1e9), o['E0'])
        nth = r['nth']
        e0 = min(x[2] for x in r['levels'] if x[0] == 0)
        epi = [x[2] for x in r['levels'] if 2 * x[0] == nth]
        if epi:
            dpi[key] = min(epi) - e0
rows = defaultdict(list)
for (model, phi, wa, L, N) in E:
    if N % 2 == 0 and (model, phi, wa, L, N + 1) in E and (model, phi, wa, L, N - 1) in E:
        d = E[(model, phi, wa, L, N + 1)] + E[(model, phi, wa, L, N - 1)] - 2 * E[(model, phi, wa, L, N)]
        rows[(model, phi, wa)].append((N / L, L, N, d))
for k in sorted(rows):
    print(f"\n{k[0]} phi={k[1]} |w|={k[2]}: parity gap at even N")
    for (n, L, N, d) in sorted(rows[k]):
        print(f"   n={n:.4f} L={L:3d} N={N:2d}: Delta_pb={d:.5f}  L*Delta_pb={L*d:.3f}")
print("\nflux-period test delta_pi = E_min(theta=pi) - E_min(theta=0) at fixed N:")
for k in sorted(dpi):
    print(f"   {k}: delta_pi={dpi[k]:+.5f}  L*delta_pi={k[3]*dpi[k]:+.3f}")
