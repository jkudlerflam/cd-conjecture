"""C186: particle-number convergence of the K = pi two-quasiparticle level at fixed (phi, mu) taken from pairbind2.py.
For N_max = 3 ... N_top: dressed one-particle levels E1 at the two valley momenta and the lowest K = pi level.
The single quasiparticle converges by N_max ~ 4 to 5 (dressing.out); the pair needs about twice the extra particles."""
import numpy as np, json, sys, time
from ring_mixed import Ring, model
src, v, Ntop = sys.argv[1], float(sys.argv[2]), int(sys.argv[3])
d = [json.loads(l) for l in open(src) if json.loads(l)['v'] == v and json.loads(l)['rmode'] == 'locked']
d = max(d, key=lambda x: x['Nmax'])
L = d['L']; m = model(d['u'], d['phi'], v, d['mu'])
ia = int(round(d['qa'] * L / (2 * np.pi))) % L; ib = int(round(d['qb'] * L / (2 * np.pi))) % L
iK = (ia + ib) % L
for Nmax in range(3, Ntop + 1):
    t0 = time.time()
    ring = Ring(L, Nmax)
    out = dict(u=d['u'], v=v, phi=d['phi'], mu=d['mu'], L=L, Nmax=Nmax)
    for name, i in (("E1a", ia), ("E1b", ib)):
        e, Nx = ring.lowest(m, 2 * np.pi * i / L, k=3)
        j = next(k for k in range(len(e)) if 0.5 < Nx[k] < 1.7)
        out[name] = float(e[j]); out[name + "_N"] = float(Nx[j])
    e, Nx = ring.lowest(m, 2 * np.pi * iK / L, k=2)
    out["Epair"] = float(e[0]); out["Npair"] = float(Nx[0])
    out["B_same_cap"] = out["E1a"] + out["E1b"] - out["Epair"]
    out["sec"] = round(time.time() - t0, 1)
    print(json.dumps(out), flush=True)
