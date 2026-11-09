"""C186: like-valley (aa, bb) and unlike (ab) two-quasiparticle channels, and the three-quasiparticle channels aab, abb,
of the mixed cusp model at the (phi, mu) tuned by pairbind2.py (read from its JSON output). Ring L, cap Nmax.
For each channel K, E_low(K) is compared with the free threshold T(K) built from the dressed one-particle levels E1."""
import numpy as np, json, sys
from ring_mixed import Ring, model
src, Nmax = sys.argv[1], int(sys.argv[2])
for line in open(src):
    d = json.loads(line)
    if d.get('rmode') != 'locked' or d['Nmax'] != Nmax:
        continue
    L = d['L']; ring = Ring(L, Nmax)
    m = model(d['u'], d['phi'], d['v'], d['mu'])
    E1 = np.full(L, np.inf)
    for i in range(L):
        e, Nx = ring.lowest(m, 2*np.pi*i/L, k=5)
        for ei, ni in zip(e, Nx):
            if 0.5 < ni < 1.7:
                E1[i] = ei; break
    ia = int(round(d['qa']*L/(2*np.pi))) % L; ib = int(round(d['qb']*L/(2*np.pi))) % L
    T2 = lambda iK: min(E1[j] + E1[(iK-j) % L] for j in range(L))
    T3 = lambda iK: min(E1[j] + E1[k] + E1[(iK-j-k) % L] for j in range(L) for k in range(L))
    out = dict(u=d['u'], v=d['v'], L=L, Nmax=Nmax, phi=d['phi'], mu=d['mu'])
    for name, iK, T in (("aa", (2*ia) % L, T2), ("bb", (2*ib) % L, T2), ("ab", (ia+ib) % L, T2),
                        ("aab", (2*ia+ib) % L, T3), ("abb", (ia+2*ib) % L, T3)):
        e, Nx = ring.lowest(m, 2*np.pi*iK/L, k=2)
        t = T(iK)
        out[name] = dict(K=2*np.pi*iK/L, E_low=float(e[0]), N=float(Nx[0]), T=float(t), B=float(t - e[0]))
    print(json.dumps(out), flush=True)
