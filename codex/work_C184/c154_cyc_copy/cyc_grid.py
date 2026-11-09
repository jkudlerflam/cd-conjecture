"""C154: scan of the fundamental triangle of the cyclic stratum, h = XZ + b YX + g ZY with 1 >= b >= |g|
(Claude Opus agent, 2026-11-07). Rings L = 12, 18 (momentum ED, C155's RingTable).
For each point: ground-state momentum, the lowest triplet (three degenerate levels at relative momenta pi + {0, +-2pi/3}),
the singlet at relative momentum pi, the current level at relative momentum 1 (x = 1 in SU(2)_1), the ratios
r_s = singlet/triplet and r_c = current/triplet (SU(2)_1 with log corrections: r_s >~ 1, r_c ~ 2), the Casimir c
from E0/L at L = 12, 18 with v from the current level at L = 18, and the half-ring entropy slope c_S."""
import json, sys, time
import numpy as np
import est_c155 as E

def bond_abg(a, b, g):
    J = np.zeros((3, 3)); J[0, 2] = a; J[1, 0] = b; J[2, 1] = g
    return E.bond_from_JB(J, np.zeros(3))

def analyse(h, L):
    d = E.ring_data(h, L, nev=3)
    lev = d['lev']; E0 = d['E0']; m0 = d['m0']
    rel = lambda dm: (m0 + dm) % L
    trip = [lev[rel(L // 2 + s), 0] if s == 0 else lev[rel(L // 2 + s), 0] for s in (0, L // 3, -L // 3)]
    pi_levels = np.sort(lev[rel(L // 2)])
    # triplet: lowest level at pi + 2pi/3; singlet: the level at pi not degenerate with it
    t = lev[rel(L // 2 + L // 3), 0]
    cand = [x for x in pi_levels if abs(x - t) > 1e-7]
    s = cand[0] if cand else np.nan
    cur = min(lev[rel(1), 0], lev[rel(-1), 0])
    S = E.entropies(d['psi'], L)
    allw = np.sort(lev.ravel())
    return dict(L=L, E0=E0, m0=int(m0), x_t=L * (t - E0), x_s=L * (s - E0), x_c=L * (cur - E0),
                trip_split=float(L * (max(trip) - min(trip))), gap1=L * (allw[1] - E0), Shalf=float(S[L // 2]))

out = sys.argv[1] if len(sys.argv) > 1 else 'cyc_grid.jsonl'
pts = []
for b in (0.1, 0.25, 0.4, 0.55, 0.7, 0.85, 1.0):
    for f in (-0.9, -0.6, -0.3, 0.0, 0.3, 0.6, 0.9):
        pts.append((1.0, b, round(f * b, 4)))
pts += [(1.0, 1.0, 1.0), (1.0, 1.0, -0.98), (1.0, 0.05, 0.0), (1.0, 0.05, 0.045), (1.0, 0.05, -0.045)]
if len(sys.argv) > 2:
    pts = [tuple(float(x) for x in p.split(',')) for p in sys.argv[2:]]
for (a, b, g) in pts:
    t0 = time.time()
    h = bond_abg(a, b, g)
    r12, r18 = analyse(h, 12), analyse(h, 18)
    cv = 6 * (r18['E0'] / 18 - r12['E0'] / 12) / (np.pi * (1 / 144 - 1 / 324))
    v = r18['x_c'] / (2 * np.pi)
    rec = dict(a=a, b=b, g=g, r12=r12, r18=r18, c_cas=cv / v if v > 1e-6 else None, v=v,
               c_S=3 * (r18['Shalf'] - r12['Shalf']) / np.log(18 / 12), sec=round(time.time() - t0, 1))
    with open(out, 'a') as fo:
        fo.write(json.dumps(rec) + '\n')
    print(f"(1,{b},{g}): m0={r12['m0']},{r18['m0']} x_t={r12['x_t']:.3f},{r18['x_t']:.3f} x_s={r12['x_s']:.3f},{r18['x_s']:.3f} "
          f"x_c={r12['x_c']:.3f},{r18['x_c']:.3f} split={r18['trip_split']:.1e} gap1(18)={r18['gap1']:.3f} "
          f"r_s={r18['x_s']/r18['x_t']:.3f} r_c={r18['x_c']/r18['x_t']:.3f} c_cas={rec['c_cas']:.3f} c_S={rec['c_S']:.3f} ({rec['sec']}s)", flush=True)
