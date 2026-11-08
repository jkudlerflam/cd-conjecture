"""C164: two-particle levels on the ring at a tuned two-valley point. For every sector m, the lowest level that is
not the one-particle state is compared with the non-interacting floor min_m1 [E1(m1) + E1(m - m1)] (lattice values).
dE2(m) < 0 at the sector where the floor is made of two A (or two B) particles signals like-valley binding;
the AA, BB and AB shifts are listed. usage: twobody.py endpoints.jsonl index L"""
import sys, json, numpy as np, est, qband
r = json.loads(open(sys.argv[1]).readlines()[int(sys.argv[2])]); L = int(sys.argv[3])
v = np.array(r['hist'][-1]['v']); inf = r['hist'][-1]['info']
J, f = est.J_f(v); h = est.bond_from_JB(J, f)
b = qband.band(h, L, nev=4); lev = b['lev'] - b['E0']; E1 = b['E1']
mA = int(round(inf['kA'] / (2 * np.pi / L))) % L; mB = int(round(inf['kB'] / (2 * np.pi / L))) % L
print(f"L={L} kA={inf['kA']:.3f} (m={mA}) kB={inf['kB']:.3f} (m={mB}); E1 =", np.round(E1, 4).tolist())
for m in range(L):
    pairs = sorted(((E1[m1] + E1[(m - m1) % L], m1, (m - m1) % L) for m1 in range(L)))
    fl, m1, m2 = pairs[0]
    ls = sorted(lev[m]); ls = [x for x in ls if x > 1e-9]
    one = E1[m]
    two = [x for x in ls if abs(x - one) > 1e-9] if one < fl else ls
    lab = ''.join(sorted(('A' if x == mA else 'B' if x == mB else '.') for x in (m1, m2)))
    print(f"m={m:2d} floor={fl:.4f} ({m1},{m2}) {lab:3s} lowest two-particle={two[0]:.4f} shift={two[0] - fl:+.4f}")
