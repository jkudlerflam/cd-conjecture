"""C154: combine ring data L = 12, 18 (cyc_spec.jsonl) and L = 24 (typhon_out/ed24.jsonl) on the YX line and at Heisenberg.
x_t, x_s, x_c = L (E - E0) of the triplet (rel. momentum pi + 2pi/3), the singlet (rel. pi, non-degenerate with the triplet)
and the current (rel. momentum 2pi/L). Velocity v = x_c / 2pi (current has x = 1 in SU(2)_1).
Casimir: E0/L = e0 - pi c v/(6 L^2); pairwise cv from (12,18) and (18,24); c = cv / v(L=24) (or v(18))."""
import json, numpy as np
spec = {}
for l in open('cyc_spec.jsonl'):
    r = json.loads(l); key = (r['a'], r['b'], r['g'])
    L = r['L']; lev = np.array(r['lev']); m0 = r['m0']
    E0 = r['E0']
    t = lev[(m0 + L // 2 + L // 3) % L, 0]; pil = np.sort(lev[(m0 + L // 2) % L]); s = [x for x in pil if abs(x - t) > 1e-6][0]
    c = min(lev[(m0 + 1) % L, 0], lev[(m0 - 1) % L, 0])
    spec.setdefault(key, {})[L] = dict(E0=E0, x_t=t, x_s=s, x_c=c)
ed = {}
for l in open('typhon_out/ed24.jsonl'):
    r = json.loads(l); key = (r['a'], r['b'], r['g']); ed.setdefault(key, {})[r['m']] = r['levels']
for key, d in sorted(ed.items()):
    if not all(m in d for m in (0, 1, 4, 12)): continue
    L = 24; E0 = d[0][0]
    t = d[4][0]; s = [x for x in d[12] if abs(x - t) > 1e-7][0]; c = d[1][0]
    spec.setdefault(key, {})[24] = dict(E0=E0, x_t=L * (t - E0), x_s=L * (s - E0), x_c=L * (c - E0))
print(" (a,b,g)      L : x_t     x_s     x_c    r_s=x_s/x_t  r_c=x_c/x_t   v=x_c/2pi")
for key in sorted(spec):
    d = spec[key]
    for L in sorted(d):
        q = d[L]
        print(f"{key} {L:3d}: {q['x_t']:7.3f} {q['x_s']:7.3f} {q['x_c']:7.3f}   {q['x_s']/q['x_t']:6.3f}   {q['x_c']/q['x_t']:6.3f}   {q['x_c']/(2*np.pi):7.4f}")
    Ls = sorted(d)
    out = []
    for L1, L2 in zip(Ls[:-1], Ls[1:]):
        cv = 6 * (d[L2]['E0'] / L2 - d[L1]['E0'] / L1) / (np.pi * (1 / L1**2 - 1 / L2**2))
        out.append(f"cv({L1},{L2})={cv:.4f} -> c={cv/(d[L2]['x_c']/(2*np.pi)):.3f}")
    if len(Ls) >= 3:
        A = np.array([[1, -np.pi / (6 * L**2)] for L in Ls[-3:]]); y = np.array([d[L]['E0'] / L for L in Ls[-3:]])
    print("     " + "; ".join(out))
